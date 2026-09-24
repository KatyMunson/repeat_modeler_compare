#!/usr/bin/env python3
"""Genome identity fingerprint (addendum §4, reworked).

write:  md5 of the prepped FASTA, sequence count, total bp, then every
        sequence's name and length (sorted by name).
check:  compare two fingerprints; exit non-zero with a readable diff if
        they differ. repeatmodeler uses this before -recoverDir-ing an
        existing RM_* directory, so a run built on a different genome
        (pre-scaffold vs scaffolded assembly, or a different
        genome_prep.test_subsample_bp) is never silently resumed.
Stdlib only."""

import argparse
import hashlib
import sys

from fasta_utils import iter_fasta, seq_id


def md5_file(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write(args):
    lengths = {}
    for header, seq in iter_fasta(args.fasta):
        lengths[seq_id(header)] = len(seq)
    with open(args.out, "w") as fh:
        fh.write(f"#md5\t{md5_file(args.fasta)}\n")
        fh.write(f"#n_seqs\t{len(lengths)}\n")
        fh.write(f"#total_bp\t{sum(lengths.values())}\n")
        for name in sorted(lengths):
            fh.write(f"{name}\t{lengths[name]}\n")


def read(path):
    meta = {}
    lengths = {}
    with open(path) as fh:
        for line in fh:
            key, value = line.rstrip("\n").split("\t")
            if key.startswith("#"):
                meta[key[1:]] = value
            else:
                lengths[key] = int(value)
    return meta, lengths


def check(args):
    meta_a, len_a = read(args.expected)
    meta_b, len_b = read(args.observed)
    if meta_a["md5"] == meta_b["md5"]:
        print(f"[fingerprint] match (md5 {meta_a['md5']})")
        return
    only_a = sorted(set(len_a) - set(len_b))
    only_b = sorted(set(len_b) - set(len_a))
    diff_len = sorted(n for n in set(len_a) & set(len_b) if len_a[n] != len_b[n])
    msg = [
        f"[fingerprint] GENOME MISMATCH: {args.expected} vs {args.observed}",
        f"  md5 {meta_a['md5']} vs {meta_b['md5']}",
        f"  n_seqs {meta_a['n_seqs']} vs {meta_b['n_seqs']}, "
        f"total_bp {meta_a['total_bp']} vs {meta_b['total_bp']}",
        f"  {len(only_a)} names only in expected (e.g. {only_a[:3]}), "
        f"{len(only_b)} only in observed (e.g. {only_b[:3]}), "
        f"{len(diff_len)} with different lengths (e.g. {diff_len[:3]})",
    ]
    print("\n".join(msg), file=sys.stderr)
    sys.exit(1)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    w = sub.add_parser("write")
    w.add_argument("--fasta", required=True)
    w.add_argument("--out", required=True)
    c = sub.add_parser("check")
    c.add_argument("--expected", required=True)
    c.add_argument("--observed", required=True)
    args = ap.parse_args()
    write(args) if args.cmd == "write" else check(args)


if __name__ == "__main__":
    main()
