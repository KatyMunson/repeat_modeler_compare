#!/usr/bin/env python3
"""Sanitize a genome FASTA for RepeatMasker/RepeatModeler: decompress if
needed, truncate/clean headers, enforce uniqueness, optionally filter short
contigs or subsample to a fixed bp budget for wiring tests. Writes the
cleaned FASTA plus a name_map.tsv (original -> sanitized) so downstream
results can be mapped back to original contig names. Stdlib only."""

import argparse
import gzip
import re
import sys

DISALLOWED = re.compile(r"[^A-Za-z0-9_.\-]")


def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def iter_fasta(fh):
    header = None
    seq_chunks = []
    for line in fh:
        line = line.rstrip("\n")
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seq_chunks)
            header = line[1:]
            seq_chunks = []
        else:
            seq_chunks.append(line.strip())
    if header is not None:
        yield header, "".join(seq_chunks)


def sanitize(header, max_len):
    # Truncate at first whitespace (drop any description).
    name = header.split()[0] if header.split() else header
    name = DISALLOWED.sub("_", name)
    if len(name) > max_len:
        name = name[:max_len]
    if not name:
        name = "seq"
    return name


def unique_name(name, used, max_len):
    if name not in used:
        used.add(name)
        return name
    n = 2
    while True:
        suffix = f"_{n}"
        if len(suffix) >= max_len:
            raise ValueError(
                f"max_header_len={max_len} too small to disambiguate "
                f"repeated sanitized name '{name}'"
            )
        candidate = name[: max_len - len(suffix)] + suffix
        if candidate not in used:
            used.add(candidate)
            return candidate
        n += 1


def write_fasta(fh, name, seq, width=60):
    fh.write(f">{name}\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i : i + width] + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True, help="input assembly FASTA (may be .gz)")
    ap.add_argument("--min-contig-len", type=int, default=0)
    ap.add_argument("--max-header-len", type=int, default=50)
    ap.add_argument("--test-subsample-bp", type=int, default=0)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-name-map", required=True)
    args = ap.parse_args()

    used_names = set()
    total_kept_bp = 0
    n_in = 0
    n_kept = 0

    with open_maybe_gzip(args.fasta) as fh_in, open(args.out_fasta, "w") as fh_out, open(
        args.out_name_map, "w"
    ) as fh_map:
        fh_map.write("original_name\tsanitized_name\n")
        for header, seq in iter_fasta(fh_in):
            n_in += 1
            if args.min_contig_len > 0 and len(seq) < args.min_contig_len:
                continue

            if args.test_subsample_bp > 0:
                remaining = args.test_subsample_bp - total_kept_bp
                if remaining <= 0:
                    break
                if len(seq) > remaining:
                    seq = seq[:remaining]

            sanitized = sanitize(header, args.max_header_len)
            sanitized = unique_name(sanitized, used_names, args.max_header_len)

            write_fasta(fh_out, sanitized, seq)
            fh_map.write(f"{header}\t{sanitized}\n")

            total_kept_bp += len(seq)
            n_kept += 1

    print(
        f"[prep_genome] read {n_in} contigs, kept {n_kept}, "
        f"{total_kept_bp} bp written to {args.out_fasta}",
        file=sys.stderr,
    )
    if n_kept == 0:
        raise ValueError(f"No contigs survived filtering for {args.fasta} — check min_contig_len/test_subsample_bp")


if __name__ == "__main__":
    main()
