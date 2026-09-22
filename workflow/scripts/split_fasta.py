#!/usr/bin/env python3
# Copied verbatim from KatyMunson/compare_assemblies_satellites
# (common/scripts/split_fasta.py, commit 27ce4ee) for repeat_compare's
# repeatmasker scatter/gather. Not modified.
"""
split_fasta.py -- snake/boustrophedon split a FASTA file's records across N
output files, by contig (record) count, not by base count.

Canonical, monorepo-shared version. History: originally written in
vendor/rhodonite/, adapted into compare_meadowlark/scripts/split_fasta.py,
then upgraded there (in trf_consensus/workflow/scripts/split_fasta.py) from
plain round-robin (record i -> output i % N, always wrapping N-1 -> 0) to a
snake/boustrophedon assignment: the target index bounces back and forth
(0,1,...,N-1,N-1,...,1,0,0,1,...) instead of always wrapping. This is that
upgraded version, now the single copy used by every stage that scatters a
FASTA for parallel RepeatMasker/TRF processing.

Still a single sequential line-scan, no sequence-length lookahead, no
pysam -- this only actually improves load balance over plain round-robin
when contigs happen to be roughly size-ordered in the source fasta already
(common for hifiasm/verkko output -- largest contigs first -- but not
guaranteed); if the input order is arbitrary, it's no worse than plain
round-robin. Not a substitute for real bp-balanced bin-packing, which
remains out of scope here.
"""
import argparse
import sys

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--infile", required=True, help="input FASTA file")
    parser.add_argument(
        "--outputs", nargs="+", required=True, help="output FASTA chunk files"
    )
    args = parser.parse_args()

    n = len(args.outputs)
    outs = [open(f, "w") for f in args.outputs]
    current = None
    next_idx = 0
    direction = 1
    with open(args.infile) as fasta:
        for line in fasta:
            if line.startswith(">"):
                current = outs[next_idx]
                next_idx += direction
                if next_idx == n:
                    next_idx = n - 1
                    direction = -1
                elif next_idx < 0:
                    next_idx = 0
                    direction = 1
            if current is None:
                sys.exit(
                    f"ERROR: {args.infile} has sequence data before its "
                    "first '>' header"
                )
            current.write(line)

    for out in outs:
        out.close()
