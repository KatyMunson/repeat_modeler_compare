#!/usr/bin/env python3
"""Discovery-round saturation check for one species' OWN-arm masking.

Groups masked bp by the RepeatModeler round that discovered each family,
parsed from the family name in the .out "matching repeat" column:
  <code>_rnd-<N>_family-<M>  -> rnd-<N>
  <code>_ltr-<N>_family-<M>  -> ltr
  anything else              -> other (Dfam export)
bp are merged per contig within each bucket. If the final round's
families still mask a meaningful share of the genome (e.g. >1%), the
RepeatModeler sampling has not saturated -- consider
repeatmodeler.extra_args: "-numAddlRounds 1" (same for every species).
Stdlib only."""

import argparse
import re

from summarize_rm import merge_intervals

RND = re.compile(r"(?:^|_)(rnd-\d+)_family-\d+")
LTR = re.compile(r"(?:^|_)ltr-\d+_family-\d+")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-file", required=True, help="own-arm RepeatMasker .out")
    ap.add_argument("--assembly-stats", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    with open(args.assembly_stats) as fh:
        header = fh.readline().strip().split("\t")
        stats = dict(zip(header, fh.readline().strip().split("\t")))
    total_bp, non_n_bp = int(stats["total_bp"]), int(stats["non_n_bp"])

    buckets = {}
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

    def order(b):
        return (0, int(b.split("-")[1])) if b.startswith("rnd-") else (1, b)

    with open(args.out, "w") as out:
        out.write("species\tbucket\tbp\tpct_total\tpct_non_n\n")
        for bucket in sorted(buckets, key=order):
            bp = sum(merge_intervals(iv) for iv in buckets[bucket].values())
            out.write(f"{args.species}\t{bucket}\t{bp}\t"
                      f"{100.0 * bp / total_bp:.4f}\t{100.0 * bp / non_n_bp:.4f}\n")


if __name__ == "__main__":
    main()
