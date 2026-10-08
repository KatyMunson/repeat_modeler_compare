#!/usr/bin/env python3
"""Run one EDTA structural caller (EDTA_raw.pl --type tir|helitron) over
one genome group in windows, each under a timeout (README "Structural
TIR / Helitron pilot"). Runs inside the EDTA container.

The group FASTA (whole scaffolds, from group_genome.py) is cut into
pieces of <= --window-size bp (long scaffolds overlap by --overlap bp),
and pieces are packed into windows of about --window-size bp. Each window
is one EDTA_raw.pl run in its own directory, with short record names
(p1, p2, ...) so EDTA renames nothing (--convert_seq_name 0). Up to
threads / --threads-per-window windows run at once. A window still
running after --timeout seconds is killed and its pieces are written to
--timeouts (contig, start, end in scaffold coordinates), the way the LTR
tools log skipped windows. A window that exits without its result file
fails the job: nothing is dropped silently.

Outputs: --out-fa, the intact candidates with headers
">{group}_{type}_{n}#{EDTA label}"; --out-tsv, one row per candidate
with its scaffold coordinates (NA when EDTA's header gives none).
Stdlib only."""

import argparse
import glob
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

from fasta_utils import iter_fasta, seq_id, write_fasta

RESULT = {"tir": "TIR.intact.raw.fa", "helitron": "Helitron.intact.raw.fa"}
# EDTA candidate names: "p3:1200..4567#DNA/DTA" (TIR-Learner) or
# "p3_1200_4567#DNA/Helitron" (HelitronScanner); the piece is the record name.
COORDS = re.compile(r"^(p\d+)[:_](\d+)(?:\.\.|_)(\d+)")


def pieces(fasta, window, overlap):
    """[(name, contig, offset, seq)]; offset is 0-based on the scaffold."""
    out = []
    for header, seq in iter_fasta(fasta):
        contig = seq_id(header)
        step = window - overlap
        start = 0
        while True:
            chunk = seq[start:start + window]
            if chunk.strip("Nn"):
                out.append((f"p{len(out) + 1}", contig, start, chunk))
            if start + window >= len(seq):
                break
            start += step
    return out


def pack(pcs, window):
    """Greedy, in order: pieces into windows of <= window bp (a piece
    longer than that is a window of its own)."""
    wins, cur, size = [], [], 0
    for p in pcs:
        if cur and size + len(p[3]) > window:
            wins.append(cur)
            cur, size = [], 0
        cur.append(p)
        size += len(p[3])
    if cur:
        wins.append(cur)
    return wins


def find_edta_raw():
    """EDTA_raw.pl: on PATH, else next to EDTA.pl, else in the share/EDTA*
    directory bioconda installs EDTA into (images may link only EDTA.pl
    into bin/)."""
    hit = shutil.which("EDTA_raw.pl")
    if hit:
        return os.path.realpath(hit)
    cands = []
    edta = shutil.which("EDTA.pl")
    if edta:
        cands.append(os.path.join(os.path.dirname(os.path.realpath(edta)), "EDTA_raw.pl"))
    for prefix in filter(None, [os.environ.get("CONDA_PREFIX"), "/usr/local", "/opt/conda"]):
        cands += sorted(glob.glob(os.path.join(prefix, "share", "EDTA*", "EDTA_raw.pl")))
    for c in cands:
        if os.path.isfile(c):
            return c
    sys.exit("[dna_te] EDTA_raw.pl not found (PATH, next to EDTA.pl, share/EDTA*): is this running in "
             "the dna_te.container image (snakemake --use-singularity)? Looked at: " + ", ".join(cands))


def run_window(i, win, args):
    name = f"w{i}"
    d = os.path.join(args.workdir, name)
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    with open(os.path.join(d, f"{name}.fa"), "w") as fh:
        for pname, _c, _o, seq in win:
            write_fasta(fh, pname, seq)
    cmd = ["timeout", "-k", "60", str(args.timeout), "perl", args.edta_raw, "--genome", f"{name}.fa",
           "--type", args.type, "--species", args.species, "--threads", str(args.threads_per_window),
           "--convert_seq_name", "0", "--overwrite", "1"]
    bp = sum(len(p[3]) for p in win)
    print(f"[dna_te] {name}: {len(win)} pieces, {bp} bp: {' '.join(cmd)}", flush=True)
    with open(os.path.join(d, "edta_raw.log"), "w") as log:
        rc = subprocess.run(cmd, cwd=d, stdout=log, stderr=subprocess.STDOUT).returncode
    if rc in (124, 137):
        print(f"[dna_te] {name}: timed out after {args.timeout} s; logged as skipped", flush=True)
        return name, None
    hits = glob.glob(os.path.join(d, "*.EDTA.raw", f"*.{RESULT[args.type]}"))
    if rc != 0 or not hits:
        tail = open(os.path.join(d, "edta_raw.log")).read()[-3000:]
        raise RuntimeError(f"{name}: EDTA_raw.pl exit {rc}, result {'found' if hits else 'missing'}\n{tail}")
    return name, list(iter_fasta(hits[0]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True, help="group FASTA (whole scaffolds)")
    ap.add_argument("--type", required=True, choices=sorted(RESULT))
    ap.add_argument("--group", required=True)
    ap.add_argument("--species", default="others", help="EDTA_raw --species (TIR-Learner model)")
    ap.add_argument("--window-size", type=int, default=5_000_000)
    ap.add_argument("--overlap", type=int, default=100_000)
    ap.add_argument("--timeout", type=int, required=True, help="seconds per window")
    ap.add_argument("--threads", type=int, required=True)
    ap.add_argument("--threads-per-window", type=int, default=4)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out-fa", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--timeouts", required=True)
    ap.add_argument("--keep-workdir", action="store_true")
    args = ap.parse_args()
    if args.overlap >= args.window_size:
        sys.exit("--overlap must be smaller than --window-size")

    args.edta_raw = find_edta_raw()
    print(f"[dna_te] EDTA_raw.pl: {args.edta_raw}")
    pcs = pieces(args.fasta, args.window_size, args.overlap)
    wins = pack(pcs, args.window_size)
    where = {p[0]: (p[1], p[2]) for p in pcs}
    n_par = max(1, args.threads // args.threads_per_window)
    print(f"[dna_te] {args.group} {args.type}: {len(pcs)} pieces in {len(wins)} windows, {n_par} at a time")
    os.makedirs(args.workdir, exist_ok=True)

    with ThreadPoolExecutor(n_par) as pool:
        results = list(pool.map(lambda iw: run_window(iw[0], iw[1], args), enumerate(wins, 1)))

    n = 0
    with open(args.out_fa, "w") as fa, open(args.out_tsv, "w") as tsv, open(args.timeouts, "w") as to:
        tsv.write("candidate\ttype\tlabel\tcontig\tstart\tend\tlength\tedta_name\n")
        to.write("tool\tcontig\tstart\tend\n")
        for (name, recs), win in zip(results, wins):
            if recs is None:
                for _p, contig, off, seq in win:
                    to.write(f"{args.type}\t{contig}\t{off}\t{off + len(seq)}\n")
                continue
            for header, seq in recs:
                n += 1
                edta = seq_id(header)
                label = edta.split("#", 1)[1] if "#" in edta else "NA"
                cid = f"{args.group}_{args.type}_{n}"
                contig, start, end = "NA", "NA", "NA"
                m = COORDS.match(edta)
                if m and m.group(1) in where:
                    contig, off = where[m.group(1)]
                    s, e = sorted((int(m.group(2)), int(m.group(3))))
                    start, end = off + s - 1, off + e
                write_fasta(fa, f"{cid}#{label}", seq)
                tsv.write(f"{cid}\t{args.type}\t{label}\t{contig}\t{start}\t{end}\t{len(seq)}\t{edta}\n")
    skipped = sum(1 for _n, r in results if r is None)
    print(f"[dna_te] {args.group} {args.type}: {n} candidates; {skipped} of {len(wins)} windows timed out")
    if not args.keep_workdir:
        shutil.rmtree(args.workdir, ignore_errors=True)


if __name__ == "__main__":
    main()
