#!/usr/bin/env python3
"""Discovery-round saturation check for one sample's OWN-arm masking.

Groups masked bp by the RepeatModeler round that discovered each family,
parsed from the family name in the .out "matching repeat" column:
  <code>_rnd-<N>_family-<M>  -> rnd-<N>
  <code>_ltr-<N>_family-<M>  -> ltr
  anything else              -> other (Dfam export)
bp are merged per contig within each bucket. If the final round's
families still mask a meaningful share of the genome (e.g. >1%), the
RepeatModeler sampling may not have saturated -- check
discovery_round_novelty.tsv (round_novelty.py), which separates new families
from refinements of earlier ones, then consider repeatmodeler.extra_rounds.

With --tandem-table (this run's family_tandem.tsv), tandem_bp is the part
of a bucket held by tandem_family families and dispersed_pct_non_n the
rest. A late round dominated by a few Mb-scale tandem arrays is not the
same signal as one still finding many dispersed TE families: only the
dispersed share says another round would add TE families. Stdlib only."""

import argparse
import re

from summarize_rm import merge_intervals

RND = re.compile(r"(?:^|_)(rnd-\d+)_family-\d+")
LTR = re.compile(r"(?:^|_)ltr-\d+_family-\d+")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-file", required=True, help="own-arm RepeatMasker .out")
    ap.add_argument("--assembly-stats", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--tandem-table", help="family_tandem.tsv of the same run: split tandem vs dispersed bp")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    with open(args.assembly_stats) as fh:
        header = fh.readline().strip().split("\t")
        stats = dict(zip(header, fh.readline().strip().split("\t")))
    total_bp, non_n_bp = int(stats["total_bp"]), int(stats["non_n_bp"])
    tandem = set()
    if args.tandem_table:
        from family_tandem import tandem_families
        tandem = tandem_families(args.tandem_table)

    buckets = {}
    tandem_buckets = {}
    with open(args.out_file) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 11 or not f[0].isdigit():
                continue
            name = f[9]
            m = RND.search(name)
            if m:
                bucket = m.group(1)
            elif LTR.search(name):
                bucket = "ltr"
            else:
                bucket = "other"
            begin, end = sorted((int(f[5]), int(f[6])))
            buckets.setdefault(bucket, {}).setdefault(f[4], []).append((begin, end))
            if name in tandem:
                tandem_buckets.setdefault(bucket, {}).setdefault(f[4], []).append((begin, end))

    def order(b):
        return (0, int(b.split("-")[1])) if b.startswith("rnd-") else (1, b)

    with open(args.out, "w") as out:
        out.write("sample\tbucket\tbp\tpct_total\tpct_non_n\ttandem_bp\tdispersed_pct_non_n\n")
        for bucket in sorted(buckets, key=order):
            bp = sum(merge_intervals(iv) for iv in buckets[bucket].values())
            if args.tandem_table:
                t_bp = sum(merge_intervals(iv) for iv in tandem_buckets.get(bucket, {}).values())
                t_cols = f"{t_bp}\t{100.0 * (bp - t_bp) / non_n_bp:.4f}"
            else:
                t_cols = "NA\tNA"
            out.write(f"{args.sample}\t{bucket}\t{bp}\t"
                      f"{100.0 * bp / total_bp:.4f}\t{100.0 * bp / non_n_bp:.4f}\t{t_cols}\n")


if __name__ == "__main__":
    main()
