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
running after --timeout seconds is killed, and a window whose EDTA_raw.pl
exits without its result file is set aside; both are written to
--timeouts (tool, contig, start, end in scaffold coordinates, reason), the
way the LTR tools log skipped windows. The reason says which: "timeout",
or "failed exit N" with the window's zlib compression ratio (ordinary
genomic DNA ~0.25-0.30; a satellite array ~0.02-0.07, e.g. a meadowlark
contig that broke TIR-Learner's TIRvish parsing) and the kept log's path.
A failed window of several pieces is first rerun one piece at a time, so
only the piece that fails again is skipped. A failed piece that compresses
below --tandem-zlib is a tandem array: expected (bird assemblies carry
large satellite scaffolds), it never stops the job. If failed pieces of
ordinary sequence add up to more than --max-failed-frac of the group's bp,
the job fails instead: that points to a tool problem, which skip rows
would otherwise hide.

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
import time
import zlib
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


def stamp():
    return time.strftime("%Y-%m-%d %H:%M:%S")


def find_edta_raw():
    """The command that runs EDTA_raw.pl, as a list. On PATH it is called by
    that name, unresolved: in the biocontainers image bin/EDTA_raw.pl links to
    a launcher that picks the EDTA script from the name it was called by, so
    resolving the link (or calling the launcher directly) runs full EDTA.pl
    instead. Off PATH: the real script next to EDTA.pl or in share/EDTA*/,
    run with perl."""
    if shutil.which("EDTA_raw.pl"):
        return ["EDTA_raw.pl"]
    cands = []
    edta = shutil.which("EDTA.pl")
    if edta:
        cands.append(os.path.join(os.path.dirname(os.path.realpath(edta)), "EDTA_raw.pl"))
    for prefix in filter(None, [os.environ.get("CONDA_PREFIX"), "/usr/local", "/opt/conda"]):
        cands += sorted(glob.glob(os.path.join(prefix, "share", "EDTA*", "EDTA_raw.pl")))
    for c in cands:
        if os.path.isfile(c) and not os.path.islink(c):
            return ["perl", c]
    sys.exit("[dna_te] EDTA_raw.pl not found (PATH, next to EDTA.pl, share/EDTA*): is this running in "
             "the dna_te.container image (snakemake --use-singularity)? Looked at: " + ", ".join(cands))


SHIM = os.path.join(os.path.dirname(os.path.abspath(__file__)), "compat", "pandas_positional")


def edta_env():
    """EDTA_raw.pl's environment: the pandas >= 3 shim first on PYTHONPATH
    (see compat/pandas_positional/sitecustomize.py)."""
    env = dict(os.environ)
    env["PYTHONPATH"] = SHIM + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    return env


def tool_versions():
    try:
        out = subprocess.run(["python3", "-c", "import sys, pandas; print(sys.version.split()[0], pandas.__version__)"],
                             capture_output=True, text=True, env=edta_env()).stdout.split()
        return f"python {out[0]}, pandas {out[1]}" if len(out) == 2 else "python/pandas: unknown"
    except OSError:
        return "python/pandas: unknown"


def run_one(name, win, args):
    """One EDTA_raw.pl run on the pieces in win. Returns (recs, reason):
    recs is None when it timed out ("timeout") or failed ("failed exit N ...")."""
    d = os.path.join(args.workdir, name)
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    with open(os.path.join(d, f"{name}.fa"), "w") as fh:
        for pname, _c, _o, seq in win:
            write_fasta(fh, pname, seq)
    cmd = ["timeout", "-k", "60", str(args.timeout)] + args.edta_raw + ["--genome", f"{name}.fa",
           "--type", args.type, "--species", args.species, "--threads", str(args.threads_per_window),
           "--convert_seq_name", "0", "--overwrite", "1"]
    bp = sum(len(p[3]) for p in win)
    t0 = time.time()
    print(f"[dna_te] {stamp()} {name} start: {len(win)} pieces, {bp} bp: {' '.join(cmd)}", flush=True)
    with open(os.path.join(d, "edta_raw.log"), "w") as log:
        rc = subprocess.run(cmd, cwd=d, stdout=log, stderr=subprocess.STDOUT, env=edta_env()).returncode
    if rc in (124, 137):
        print(f"[dna_te] {stamp()} {name}: timed out after {args.timeout} s", flush=True)
        return None, "timeout"
    hits = glob.glob(os.path.join(d, "*.EDTA.raw", f"*.{RESULT[args.type]}"))
    if rc != 0 or not hits:
        log_path = os.path.join(d, "edta_raw.log")
        with open(log_path, errors="replace") as fh:
            tail = "".join(fh.readlines()[-30:])
        ratio = compress_ratio(win)
        kind = "tandem array" if ratio < args.tandem_zlib else "ordinary sequence"
        reason = f"failed exit {rc} ({kind}, zlib {ratio:.2f}; log {log_path})"
        print(f"[dna_te] {stamp()} {name}: {reason}\n--- last 30 lines ---\n{tail}", flush=True)
        return None, reason
    recs = list(iter_fasta(hits[0]))
    print(f"[dna_te] {stamp()} {name} done in {(time.time() - t0) / 60:.1f} min: {len(recs)} candidates", flush=True)
    shutil.rmtree(d, ignore_errors=True)
    return recs, None


def run_window(i, win, args):
    """[(pieces, recs, reason)] for window i. A window of several pieces that
    fails (not a timeout) is rerun one piece at a time, so one bad scaffold
    (e.g. a satellite array) costs only itself, not its window-mates."""
    name = f"w{i}"
    recs, reason = run_one(name, win, args)
    if recs is not None or reason == "timeout" or len(win) == 1:
        return [(win, recs, reason)]
    print(f"[dna_te] {stamp()} {name}: rerunning its {len(win)} pieces one at a time", flush=True)
    return [(([p],) + run_one(f"{name}_{p[0]}", [p], args)) for p in win]


def compress_ratio(win):
    """zlib ratio of a window's sequence (sampled: up to 2 Mb from its start)."""
    seq = "".join(p[3] for p in win)[:2_000_000].upper().encode()
    return len(zlib.compress(seq, 6)) / len(seq) if seq else 1.0


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
    ap.add_argument("--max-failed-frac", type=float, default=0.05,
                    help="fail the job if more than this fraction of the group's bp failed on ordinary "
                         "sequence (timeouts and tandem arrays don't count)")
    ap.add_argument("--tandem-zlib", type=float, default=0.1,
                    help="a failed piece compressing below this zlib ratio is a tandem array (expected, "
                         "doesn't count toward --max-failed-frac)")
    ap.add_argument("--keep-workdir", action="store_true")
    args = ap.parse_args()
    if args.overlap >= args.window_size:
        sys.exit("--overlap must be smaller than --window-size")

    args.edta_raw = find_edta_raw()
    print(f"[dna_te] EDTA_raw.pl: {' '.join(args.edta_raw)} ({shutil.which(args.edta_raw[-1]) or args.edta_raw[-1]}); {tool_versions()} (pandas >= 3: positional shim on)")
    pcs = pieces(args.fasta, args.window_size, args.overlap)
    wins = pack(pcs, args.window_size)
    where = {p[0]: (p[1], p[2]) for p in pcs}
    n_par = max(1, args.threads // args.threads_per_window)
    print(f"[dna_te] {args.group} {args.type}: {len(pcs)} pieces in {len(wins)} windows, {n_par} at a time")
    os.makedirs(args.workdir, exist_ok=True)

    with ThreadPoolExecutor(n_par) as pool:
        results = [o for outs in pool.map(lambda iw: run_window(iw[0], iw[1], args), enumerate(wins, 1))
                   for o in outs]

    n = 0
    with open(args.out_fa, "w") as fa, open(args.out_tsv, "w") as tsv, open(args.timeouts, "w") as to:
        tsv.write("candidate\ttype\tlabel\tcontig\tstart\tend\tlength\tedta_name\n")
        to.write("tool\tcontig\tstart\tend\treason\n")
        for win, recs, reason in results:
            if recs is None:
                for _p, contig, off, seq in win:
                    to.write(f"{args.type}\t{contig}\t{off}\t{off + len(seq)}\t{reason}\n")
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
    total_bp = sum(len(p[3]) for p in pcs)
    timed_out = sum(len(p[3]) for win, _r, why in results if why == "timeout" for p in win)
    failed = [(win, why) for win, _r, why in results if why and why != "timeout"]
    tandem_bp = sum(len(p[3]) for win, why in failed if "(tandem array," in why for p in win)
    ordinary = [(win, why) for win, why in failed if "(tandem array," not in why]
    ordinary_bp = sum(len(p[3]) for win, _w in ordinary for p in win)
    print(f"[dna_te] {args.group} {args.type}: {n} candidates; skipped {timed_out} bp (timeouts), "
          f"{tandem_bp} bp (failed, tandem arrays), {ordinary_bp} bp (failed, ordinary sequence), of {total_bp} bp")
    if total_bp and ordinary_bp / total_bp > args.max_failed_frac:
        sys.exit(f"[dna_te] {ordinary_bp / total_bp:.1%} of the group's bp failed on ordinary (not tandem-array) "
                 f"sequence (> --max-failed-frac {args.max_failed_frac}): likely a tool problem, not the "
                 f"sequence; failing the job. Logs kept under {args.workdir}:\n"
                 + "\n".join(f"  {','.join(p[1] for p in win)}: {why}" for win, why in ordinary))
    # successful runs removed their own directories; failed ones keep their logs
    if not args.keep_workdir and not failed:
        shutil.rmtree(args.workdir, ignore_errors=True)


if __name__ == "__main__":
    main()
