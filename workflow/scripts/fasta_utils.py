#!/usr/bin/env python3
"""Tiny shared FASTA helpers for the satellite / LTR side-pipeline scripts.
Stdlib only. (The original scripts each keep their own iter_fasta; new
scripts import this instead of adding more copies. Python puts a script's
own directory on sys.path, so `from fasta_utils import ...` works when a
script is run as `python3 workflow/scripts/<name>.py`.)"""


def iter_fasta(path):
    """Yield (header_without_gt, sequence) for every record."""
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


def seq_id(header):
    """First whitespace-delimited token of a FASTA header."""
    return header.split()[0] if header.split() else header


def write_fasta(fh, name, seq, width=60):
    fh.write(f">{name}\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i : i + width] + "\n")


def n_runs(seq):
    """0-based half-open [start, end) runs of N/n in seq."""
    runs = []
    i = 0
    n = len(seq)
    upper = seq.upper()
    while i < n:
        j = upper.find("N", i)
        if j < 0:
            break
        k = j
        while k < n and upper[k] == "N":
            k += 1
        runs.append((j, k))
        i = k
    return runs


def merge_half_open(intervals):
    """Merge 0-based half-open intervals; returns a sorted, disjoint list."""
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            if end > merged[-1][1]:
                merged[-1] = (merged[-1][0], end)
        else:
            merged.append((start, end))
    return merged


def intersect_bp(a, b):
    """bp shared by two sorted, disjoint half-open interval lists."""
    i = j = 0
    total = 0
    while i < len(a) and j < len(b):
        lo = max(a[i][0], b[j][0])
        hi = min(a[i][1], b[j][1])
        if lo < hi:
            total += hi - lo
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return total
