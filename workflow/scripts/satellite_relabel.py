#!/usr/bin/env python3
"""Relabel de novo "Unknown" families that are really curated satellites
(instead of clustering the satellite library against them with cd-hit,
which has no notion of rotation/phase and could discard the informative
harmonized name).

build: write
   - targets.fa: every satellite-library monomer concatenated x3, so a
     consensus that starts at a different phase/rotation of the repeat
     still aligns end to end (RepeatMasker then aligns the family
     consensi against these tandem targets);
   - queries.fa: the classified families renamed q1..qN (RepeatMasker
     limits sequence names to 50 characters and dislikes '#'), plus a
     queries.tsv map back to the original headers.
apply: parse RepeatMasker's .out for queries.fa and, for each family
   classified "Unknown" whose consensus is covered >= --min-cov by
   satellite hits, rewrite its label to "#Satellite". The best-covering
   satellite name is recorded in the TSV, never written into the family
   name, so harmonized names are never renamed or dropped. Families
   with another class and high satellite coverage are reported only.
Stdlib only."""

import argparse
import sys

from fasta_utils import iter_fasta, write_fasta


def build(args):
    n = 0
    with open(args.targets, "w") as out:
        for header, seq in iter_fasta(args.satellite_lib):
            name = header.split()[0].split("#", 1)[0]
            cls = header.split()[0].split("#", 1)[1] if "#" in header.split()[0] else "Satellite"
            write_fasta(out, f"{name}#{cls}", seq * 3)
            n += 1
    if n == 0:
        raise ValueError(f"no sequences in {args.satellite_lib}")
    with open(args.queries, "w") as out, open(args.map, "w") as mp:
        mp.write("query\theader\tlength\n")
        for i, (header, seq) in enumerate(iter_fasta(args.families), 1):
            write_fasta(out, f"q{i}", seq)
            mp.write(f"q{i}\t{header}\t{len(seq)}\n")
    print(f"[satellite_relabel] {n} satellite targets (x3)", file=sys.stderr)


def parse_out(path):
    """{query: [(begin, end, repeat_name)]} from a RepeatMasker .out."""
    hits = {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 11 or not f[0].isdigit():
                continue
            begin, end = sorted((int(f[5]), int(f[6])))
            hits.setdefault(f[4], []).append((begin, end, f[9]))
    return hits


def merged_bp(intervals):
    total = 0
    cur_s = cur_e = None
    for s, e in sorted(intervals):
        if cur_e is not None and s <= cur_e + 1:
            cur_e = max(cur_e, e)
            continue
        if cur_e is not None:
            total += cur_e - cur_s + 1
        cur_s, cur_e = s, e
    if cur_e is not None:
        total += cur_e - cur_s + 1
    return total


def apply(args):
    headers = {}
    with open(args.map) as fh:
        fh.readline()
        for line in fh:
            q, header, length = line.rstrip("\n").split("\t")
            headers[q] = (header, int(length))
    hits = parse_out(args.out_file)

    decisions = {}
    with open(args.tsv, "w") as tsv:
        tsv.write("family\tclass_family\tlength\tsatellite_bp\tfrac\tbest_satellite\taction\n")
        for q, (header, length) in headers.items():
            token = header.split()[0]
            name, cls = token.split("#", 1) if "#" in token else (token, "Unknown")
            qh = hits.get(q, [])
            if not qh:
                continue
            frac = merged_bp([(b, e) for b, e, _ in qh]) / length if length else 0.0
            by_name = {}
            for b, e, rep in qh:
                by_name.setdefault(rep, []).append((b, e))
            best = max(by_name, key=lambda r: merged_bp(by_name[r]))
            if frac >= args.min_cov and cls == "Unknown":
                action = "relabeled"
                decisions[q] = f"{name}#Satellite"
            elif frac >= args.min_cov:
                action = "kept_class"
            else:
                action = "below_min_cov"
            tsv.write(f"{name}\t{cls}\t{length}\t{merged_bp([(b, e) for b, e, _ in qh])}\t{frac:.3f}\t{best}\t{action}\n")

    n = 0
    with open(args.out_fasta, "w") as out:
        for i, (header, seq) in enumerate(iter_fasta(args.families), 1):
            q = f"q{i}"
            if headers[q][0] != header:
                raise ValueError(f"{args.families} changed since build ({q})")
            if q in decisions:
                parts = header.split(None, 1)
                header = decisions[q] + (" " + parts[1] if len(parts) > 1 else "")
                n += 1
            write_fasta(out, header, seq)
    print(f"[satellite_relabel] relabeled {n} Unknown families as Satellite (min_cov {args.min_cov})",
          file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("--satellite-lib", required=True)
    b.add_argument("--families", required=True)
    b.add_argument("--targets", required=True)
    b.add_argument("--queries", required=True)
    b.add_argument("--map", required=True)
    a = sub.add_parser("apply")
    a.add_argument("--families", required=True)
    a.add_argument("--map", required=True)
    a.add_argument("--out-file", required=True, help="RepeatMasker .out for queries.fa")
    a.add_argument("--min-cov", type=float, default=0.5)
    a.add_argument("--out-fasta", required=True)
    a.add_argument("--tsv", required=True)
    args = ap.parse_args()
    build(args) if args.cmd == "build" else apply(args)


if __name__ == "__main__":
    main()
