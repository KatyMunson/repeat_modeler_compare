#!/usr/bin/env python3
"""Hard-mask (N) BED intervals in a FASTA, for LTR discovery only (the
masked genome is never used for anything else). BED is 0-based half-open.
Sequences with no BED intervals are copied unchanged. Stdlib only, so no
bedtools module is needed (see README "Environment / module policy")."""

import argparse
import sys

from fasta_utils import iter_fasta, merge_half_open, seq_id, write_fasta


def read_bed(path):
    by_contig = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            by_contig.setdefault(fields[0], []).append((int(fields[1]), int(fields[2])))
    return {c: merge_half_open(iv) for c, iv in by_contig.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--bed", help="intervals to mask; omit to copy the FASTA unmasked")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    intervals = read_bed(args.bed) if args.bed else {}
    seen = set()
    masked_bp = 0
    total_bp = 0
    with open(args.out, "w") as out:
        for header, seq in iter_fasta(args.fasta):
            name = seq_id(header)
            seen.add(name)
            total_bp += len(seq)
            ivs = intervals.get(name)
            if ivs:
                buf = bytearray(seq, "ascii")
                for start, end in ivs:
                    if end > len(buf):
                        raise ValueError(f"BED interval {name}:{start}-{end} beyond sequence length {len(buf)}")
                    buf[start:end] = b"N" * (end - start)
                    masked_bp += end - start
                seq = buf.decode("ascii")
            write_fasta(out, name, seq)
    missing = sorted(set(intervals) - seen)
    if missing:
        raise ValueError(f"BED names not in FASTA (wrong genome?): {missing[:5]} ({len(missing)} total)")
    pct = 100.0 * masked_bp / total_bp if total_bp else 0.0
    print(f"[mask_bed] masked {masked_bp} of {total_bp} bp ({pct:.2f}%)", file=sys.stderr)


if __name__ == "__main__":
    main()
