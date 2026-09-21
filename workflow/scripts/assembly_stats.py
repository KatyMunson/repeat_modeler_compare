#!/usr/bin/env python3
"""Compute basic assembly-quality covariates (contig count, total/N/non-N bp,
N50, largest contig) for one species' prepped genome FASTA. These are the
denominators used everywhere downstream (class_composition.tsv etc). Stdlib
only."""

import argparse


def iter_fasta(path):
    header = None
    seq_chunks = []
    with open(path) as fh:
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


def n50(lengths):
    lengths = sorted(lengths, reverse=True)
    total = sum(lengths)
    if total == 0:
        return 0
    half = total / 2.0
    cumulative = 0
    for length in lengths:
        cumulative += length
        if cumulative >= half:
            return length
    return lengths[-1] if lengths else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True, help="prepped genome FASTA (already sanitized)")
    ap.add_argument("--species-id", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    lengths = []
    total_bp = 0
    n_bp = 0

    for _header, seq in iter_fasta(args.fasta):
        length = len(seq)
        lengths.append(length)
        total_bp += length
        n_bp += seq.upper().count("N")

    non_n_bp = total_bp - n_bp
    contig_count = len(lengths)
    largest = max(lengths) if lengths else 0
    contig_n50 = n50(lengths)

    with open(args.out, "w") as fh:
        fh.write(
            "species_id\tcontig_count\ttotal_bp\tn_bp\tnon_n_bp\tcontig_n50\tlargest_contig\n"
        )
        fh.write(
            f"{args.species_id}\t{contig_count}\t{total_bp}\t{n_bp}\t{non_n_bp}\t"
            f"{contig_n50}\t{largest}\n"
        )


if __name__ == "__main__":
    main()
