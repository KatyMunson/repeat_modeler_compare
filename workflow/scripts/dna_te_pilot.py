#!/usr/bin/env python3
"""Structural TIR / Helitron pilot (README "Structural TIR / Helitron
pilot"). Report-only: the library is not changed.

  gather   per-group candidates from dna_te_windows.py -> one FASTA per type
           (candidates reported twice by overlapping windows dropped),
           candidates.tsv, and skipped_windows.tsv (timed-out windows,
           merged per tool, as ltr/skipped_windows.tsv).
  pilot    per sample and type: cd-hit-est the candidates (library.cdhit
           settings, both strands); cd-hit-est-2d the cluster
           representatives against the sample's own -families.fa and
           against the shared library; blastn (dc-megablast) the
           unmatched representatives against the genome. Writes one row
           per type and a per-cluster table.
  combine  per-sample rows -> summary/dna_te_pilot.tsv, with the decision
           rule on its first line.

The bp the unmatched clusters would mask is a blastn proxy: the union of
their dc-megablast HSPs on the genome, which is less sensitive than
RepeatMasker -s (old, diverged copies are missed), so it is a lower bound.
`new_bp` is the part of it the shared-arm .out does not mask already; the
decision uses new_bp. Tools (cd-hit-est, cd-hit-est-2d, makeblastdb,
blastn) come from PATH. Stdlib otherwise."""

import argparse
import os
import re
import subprocess
import sys

from fasta_utils import blast_subject_resolver, iter_fasta, merge_half_open, seq_id, write_fasta
from summarize_rm import read_assembly_stats

TYPES = ("tir", "helitron")
# EDTA superfamily codes -> the RepeatMasker / Dfam names in our libraries
SUPERFAMILY = {"DTA": "hAT", "DTC": "CMC", "DTH": "PIF-Harbinger", "DTM": "MULE", "DTT": "TcMar",
               "Helitron": "Helitron"}
SF_CODE = re.compile(r"(DT[ACHMT]|Helitron)")
BLAST_FMT = "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"


# ---------------------------------------------------------------- gather


def gather(args):
    os.makedirs(args.outdir, exist_ok=True)
    seen = set()
    n_in = {t: 0 for t in TYPES}
    n_out = {t: 0 for t in TYPES}
    fas = {t: open(os.path.join(args.outdir, f"{t}.candidates.fa"), "w") for t in TYPES}
    with open(os.path.join(args.outdir, "candidates.tsv"), "w") as out:
        header = None
        for tsv in args.tsv:
            fa = tsv[: -len(".tsv")] + ".fa"
            seqs = {seq_id(h).split("#", 1)[0]: (seq_id(h), s) for h, s in iter_fasta(fa)}
            with open(tsv) as fh:
                h = fh.readline()
                if header is None:
                    header = h
                    out.write(h)
                cols = h.rstrip("\n").split("\t")
                for line in fh:
                    r = dict(zip(cols, line.rstrip("\n").split("\t")))
                    t = r["type"]
                    n_in[t] += 1
                    key = (t, r["contig"], r["start"], r["end"])
                    if r["contig"] != "NA" and key in seen:
                        continue  # same element from two overlapping windows
                    seen.add(key)
                    name, seq = seqs[r["candidate"]]
                    write_fasta(fas[t], f"{args.sample}_{name}", seq)
                    out.write(line.replace(r["candidate"], f"{args.sample}_{r['candidate']}", 1))
                    n_out[t] += 1
    for fh in fas.values():
        fh.close()
    for t in TYPES:
        print(f"[dna_te] {args.sample} {t}: {n_in[t]} candidates, {n_out[t]} after dropping window-overlap duplicates")

    skipped = {}
    for path in args.timeouts:
        with open(path) as fh:
            fh.readline()
            for line in fh:
                tool, contig, s, e = line.rstrip("\n").split("\t")
                skipped.setdefault((tool, contig), []).append((int(s), int(e)))
    with open(os.path.join(args.outdir, "skipped_windows.tsv"), "w") as out:
        out.write("tool\tcontig\tstart\tend\tbp\n")
        for (tool, contig), ivs in sorted(skipped.items()):
            for s, e in merge_half_open(ivs):
                out.write(f"{tool}\t{contig}\t{s}\t{e}\t{e - s}\n")
                print(f"[dna_te] skipped ({tool} timed out): {contig}:{s}-{e}")


# ---------------------------------------------------------------- pilot


def run(cmd):
    print("[dna_te] " + " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def read_clstr(path):
    """[[(name, is_rep)]] from a cd-hit .clstr; names without '>' and '...'."""
    clusters = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">Cluster"):
                clusters.append([])
                continue
            m = re.search(r">(.+?)\.\.\.", line)
            if m:
                clusters[-1].append((m.group(1), line.rstrip().endswith("*")))
    return clusters


def label_of(name):
    return name.split("#", 1)[1] if "#" in name else "NA"


def lib_names(path):
    return {seq_id(h) for h, _s in iter_fasta(path)}


def cdhit_args(args):
    return ["-c", str(args.identity), "-aS", str(args.coverage_short), "-n", str(args.word_size),
            "-G", "0", "-d", "0", "-r", "1", "-M", "0", "-T", str(args.threads)]


def match_library(args, reps_fa, lib_fa, tag, wd):
    """{rep: library family it clusters with} via cd-hit-est-2d."""
    out = os.path.join(wd, f"vs_{tag}")
    run(["cd-hit-est-2d", "-i", lib_fa, "-i2", reps_fa, "-o", out] + cdhit_args(args))
    lib = lib_names(lib_fa)
    matched = {}
    for cl in read_clstr(out + ".clstr"):
        fams = [n for n, _r in cl if n in lib]
        if fams:
            for n, _r in cl:
                if n not in lib:
                    matched[n] = fams[0]
    return matched


def blast_cover(args, reps_fa, wd, tag):
    """{contig: merged [(s, e)]} covered by HSPs of reps_fa on the genome."""
    out = os.path.join(wd, f"{tag}.unmatched_vs_genome.tsv")
    run(["blastn", "-task", "dc-megablast", "-query", reps_fa, "-db", os.path.join(wd, "genome"),
         "-evalue", str(args.evalue), "-max_target_seqs", str(args.max_targets), "-num_threads", str(args.threads),
         "-outfmt", BLAST_FMT, "-out", out])
    ivs = {}
    with open(out) as fh:
        for line in fh:
            f = line.split("\t")
            s, e = sorted((int(f[6]), int(f[7])))
            ivs.setdefault(args.resolve(f[1]), []).append((s - 1, e))
    return {c: merge_half_open(v) for c, v in ivs.items()}


def out_masked(path):
    ivs = {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 11 or not f[0].isdigit():
                continue
            try:
                b, e = sorted((int(f[5]), int(f[6])))
            except ValueError:
                continue
            ivs.setdefault(f[4], []).append((b - 1, e))
    return {c: merge_half_open(v) for c, v in ivs.items()}


def uncovered_bp(cover, masked):
    """bp of cover outside masked (both {contig: sorted disjoint})."""
    total = 0
    for c, ivs in cover.items():
        m = masked.get(c, [])
        j = 0
        for s, e in ivs:
            bp = e - s
            while j < len(m) and m[j][1] <= s:
                j += 1
            k = j
            while k < len(m) and m[k][0] < e:
                bp -= max(0, min(e, m[k][1]) - max(s, m[k][0]))
                k += 1
            total += bp
    return total


def superfamily(label):
    m = SF_CODE.search(label)
    return m.group(1) if m else None


def pilot(args):
    wd = args.workdir
    os.makedirs(wd, exist_ok=True)
    with open(args.fingerprint) as fh:
        args.resolve = blast_subject_resolver(l.split("\t", 1)[0] for l in fh if not l.startswith("#"))
    _total, non_n = read_assembly_stats(args.assembly_stats)
    shared_classes = {label_of(n) for n in lib_names(args.shared_library)}
    masked = out_masked(args.out_file)
    genome_db = False

    rows = []
    with open(args.clusters_out, "w") as cl_out:
        cl_out.write("sample\ttype\tcluster_rep\tlabel\tn_members\town_match\town_class_family\t"
                     "shared_match\tshared_class_family\n")
        for t, cand in zip(TYPES, (args.tir, args.helitron)):
            n_cand = sum(1 for _ in iter_fasta(cand))
            row = {"sample": args.sample, "type": t, "candidates": n_cand}
            if not n_cand:
                row.update(clusters=0, matched_own=0, matched_shared=0, matched=0, unmatched=0,
                           unmatched_cover_bp=0, unmatched_new_bp=0, new_pct_non_n="0.0000",
                           matched_classes=".", absent_superfamilies=".")
                rows.append(row)
                continue
            nr = os.path.join(wd, f"{t}.nr.fa")
            run(["cd-hit-est", "-i", cand, "-o", nr] + cdhit_args(args))
            clusters = read_clstr(nr + ".clstr")
            reps = {n: len(cl) for cl in clusters for n, r in cl if r}
            own = match_library(args, nr, args.own_families, f"{t}_own", wd)
            shared = match_library(args, nr, args.shared_library, f"{t}_shared", wd)
            unmatched = [r for r in reps if r not in own and r not in shared]
            classes = {}
            for r in reps:
                fam = shared.get(r) or own.get(r)
                if fam:
                    cf = label_of(fam)
                    classes[cf] = classes.get(cf, 0) + 1
                cl_out.write(f"{args.sample}\t{t}\t{r}\t{label_of(r)}\t{reps[r]}\t{own.get(r, '.')}\t"
                             f"{label_of(own[r]) if r in own else '.'}\t{shared.get(r, '.')}\t"
                             f"{label_of(shared[r]) if r in shared else '.'}\n")
            absent = sorted({sf for r in unmatched for sf in [superfamily(label_of(r))]
                             if sf and not any(SUPERFAMILY[sf] in c for c in shared_classes)})
            cover_bp = new_bp = 0
            if unmatched:
                um = os.path.join(wd, f"{t}.unmatched.fa")
                keep = set(unmatched)
                with open(um, "w") as fh:
                    for h, s in iter_fasta(nr):
                        if seq_id(h) in keep:
                            write_fasta(fh, seq_id(h), s)
                if not genome_db:
                    run(["makeblastdb", "-in", args.genome, "-dbtype", "nucl", "-parse_seqids",
                         "-out", os.path.join(wd, "genome")])
                    genome_db = True
                cover = blast_cover(args, um, wd, t)
                cover_bp = sum(e - s for ivs in cover.values() for s, e in ivs)
                new_bp = uncovered_bp(cover, masked)
            row.update(clusters=len(reps), matched_own=sum(1 for r in reps if r in own),
                       matched_shared=sum(1 for r in reps if r in shared),
                       matched=len(reps) - len(unmatched), unmatched=len(unmatched),
                       unmatched_cover_bp=cover_bp, unmatched_new_bp=new_bp,
                       new_pct_non_n=f"{100.0 * new_bp / non_n:.4f}" if non_n else "NA",
                       matched_classes=",".join(f"{c}:{n}" for c, n in sorted(classes.items(), key=lambda x: -x[1])) or ".",
                       absent_superfamilies=",".join(SUPERFAMILY[s] for s in absent) or ".")
            print(f"[dna_te] {args.sample} {t}: {row}")
            rows.append(row)

    with open(args.out, "w") as out:
        out.write("\t".join(PILOT_COLS) + "\n")
        for r in rows:
            out.write("\t".join(str(r[c]) for c in PILOT_COLS) + "\n")


PILOT_COLS = ["sample", "type", "candidates", "clusters", "matched_own", "matched_shared", "matched",
              "unmatched", "unmatched_cover_bp", "unmatched_new_bp", "new_pct_non_n", "matched_classes",
              "absent_superfamilies"]


# ---------------------------------------------------------------- combine


def combine(args):
    with open(args.out, "w") as out:
        out.write(f"# decision rule: integrate (Stage B2) if the unmatched clusters add "
                  f">= {args.min_new_pct_non_n} % of non-N bp the shared library does not mask yet "
                  f"(new_pct_non_n), or carry a superfamily absent from the shared library "
                  f"(absent_superfamilies). new bp is a blastn lower bound.\n")
        out.write("\t".join(PILOT_COLS + ["decision"]) + "\n")
        for path in args.chunks:
            with open(path) as fh:
                cols = fh.readline().rstrip("\n").split("\t")
                for line in fh:
                    r = dict(zip(cols, line.rstrip("\n").split("\t")))
                    pct = float(r["new_pct_non_n"]) if r["new_pct_non_n"] != "NA" else 0.0
                    why = []
                    if pct >= args.min_new_pct_non_n:
                        why.append(f"adds {pct:.2f}%")
                    if r["absent_superfamilies"] != ".":
                        why.append(f"new superfamily {r['absent_superfamilies']}")
                    dec = "integrate: " + "; ".join(why) if why else "skip"
                    out.write("\t".join(r[c] for c in PILOT_COLS) + f"\t{dec}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    g = sub.add_parser("gather")
    g.add_argument("--sample", required=True)
    g.add_argument("--tsv", nargs="+", required=True, help="dna_te_windows.py --out-tsv files (.fa next to each)")
    g.add_argument("--timeouts", nargs="*", default=[])
    g.add_argument("--outdir", required=True)

    p = sub.add_parser("pilot")
    p.add_argument("--sample", required=True)
    p.add_argument("--tir", required=True)
    p.add_argument("--helitron", required=True)
    p.add_argument("--own-families", required=True)
    p.add_argument("--shared-library", required=True)
    p.add_argument("--genome", required=True)
    p.add_argument("--fingerprint", required=True, help="genome fingerprint.tsv (sequence names)")
    p.add_argument("--out-file", required=True, help="shared-arm .out (bp already masked)")
    p.add_argument("--assembly-stats", required=True)
    p.add_argument("--identity", type=float, default=0.8)
    p.add_argument("--coverage-short", type=float, default=0.8)
    p.add_argument("--word-size", type=int, default=5)
    p.add_argument("--evalue", type=float, default=1e-10)
    p.add_argument("--max-targets", type=int, required=True)
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--workdir", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--clusters-out", required=True)

    c = sub.add_parser("combine")
    c.add_argument("--chunks", nargs="+", required=True)
    c.add_argument("--min-new-pct-non-n", type=float, default=0.5)
    c.add_argument("--out", required=True)

    args = ap.parse_args()
    {"gather": gather, "pilot": pilot, "combine": combine}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
