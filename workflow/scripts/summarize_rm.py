#!/usr/bin/env python3
"""Parse one (arm, species) RepeatMasker run into per-class and per-Class/Family
non-overlapping bp tables, cross-checked against the .tbl summary, plus a
long-format divergence-landscape chunk from the .divsum file. Stdlib only,
no pandas. Output rows already carry arm/species/tissue so the per-run
chunks this writes can be concatenated as-is by combine_summaries.py.

Non-overlapping bp: RepeatMasker's .out can contain multiple overlapping
hits assigned to the same class (fragmented alignments of one element, or a
lower-scoring competing call). Within one collapsed class (or one raw
Class/Family, for the finer table), intervals are merged per contig before
summing, so a base already counted for that class isn't counted again —
this is a sweep over merged intervals, not a global "each base counted once
across all classes" merge (a base can still count toward two different
classes if RepeatMasker reported an overlap between them; the .tbl
cross-check below uses the true global merge to sanity-check total masked
bp, independent of any per-class double-count).
"""

import argparse
import sys

CANONICAL_CLASSES = {
    "DNA",
    "LINE",
    "SINE",
    "LTR",
    "RC",
    "Satellite",
    "Simple_repeat",
    "Low_complexity",
    "Retroposon",
    "Unknown",
}


def collapse_class(raw_class_family):
    class_part = raw_class_family.split("/", 1)[0]
    if class_part in CANONICAL_CLASSES:
        return class_part
    print(
        f"[summarize_rm] WARNING: unmapped RepeatMasker class '{raw_class_family}' -> Other",
        file=sys.stderr,
    )
    return "Other"


def parse_out_file(path):
    """Yield (query_seq, begin, end, class_family) for every hit line."""
    with open(path) as fh:
        for line in fh:
            fields = line.split()
            if len(fields) < 11:
                continue
            if not fields[0].isdigit():
                continue  # header/blank line
            query_seq = fields[4]
            try:
                begin = int(fields[5])
                end = int(fields[6])
            except ValueError:
                continue
            class_family = fields[10]
            if begin > end:
                begin, end = end, begin
            yield query_seq, begin, end, class_family


def merge_intervals(intervals):
    """intervals: list of (begin, end), 1-based inclusive. Returns merged bp."""
    if not intervals:
        return 0
    intervals = sorted(intervals)
    bp = 0
    cur_start, cur_end = intervals[0]
    for start, end in intervals[1:]:
        if start <= cur_end + 1:
            cur_end = max(cur_end, end)
        else:
            bp += cur_end - cur_start + 1
            cur_start, cur_end = start, end
    bp += cur_end - cur_start + 1
    return bp


def bp_per_group(hits, key_fn):
    """hits: list of (query_seq, begin, end, class_family). key_fn maps a hit
    to the class bucket. Returns {bucket: bp} merged per (query_seq, bucket)."""
    by_contig_bucket = {}
    for query_seq, begin, end, class_family in hits:
        bucket = key_fn(class_family)
        by_contig_bucket.setdefault((query_seq, bucket), []).append((begin, end))

    bp_by_bucket = {}
    for (_query_seq, bucket), intervals in by_contig_bucket.items():
        bp_by_bucket[bucket] = bp_by_bucket.get(bucket, 0) + merge_intervals(intervals)
    return bp_by_bucket


def total_masked_bp(hits):
    by_contig = {}
    for query_seq, begin, end, _class_family in hits:
        by_contig.setdefault(query_seq, []).append((begin, end))
    return sum(merge_intervals(intervals) for intervals in by_contig.values())


def parse_tbl_masked_bp(path):
    with open(path) as fh:
        for line in fh:
            if "bases masked" in line:
                # e.g. "bases masked:   12345678 bp ( 10.00 %)"
                for tok in line.replace(":", " ").split():
                    if tok.isdigit():
                        return int(tok)
    return None


def read_assembly_stats(path):
    with open(path) as fh:
        header = fh.readline().strip().split("\t")
        row = fh.readline().strip().split("\t")
    d = dict(zip(header, row))
    return int(d["total_bp"]), int(d["non_n_bp"])


def parse_divsum(path, landscape_max_div):
    """Parse the 'Coverage for each repeat class and divergence (Kimura)'
    table from a RepeatMasker calcDivergenceFromAlign.pl .divsum file.
    Returns list of (kimura_bin, collapsed_class, bp). Defensive: this
    section's exact column layout has varied across RepeatMasker versions,
    so unparseable rows are skipped with a warning rather than crashing —
    verify against a real .divsum during the wiring test (spec §10.2)."""
    rows = []
    try:
        with open(path) as fh:
            lines = fh.readlines()
    except FileNotFoundError:
        print(f"[summarize_rm] WARNING: no divsum file at {path}, skipping landscape", file=sys.stderr)
        return rows

    start = None
    for i, line in enumerate(lines):
        if "Coverage for each repeat class and divergence" in line:
            start = i + 1
            break
    if start is None:
        print(f"[summarize_rm] WARNING: no Kimura coverage table found in {path}", file=sys.stderr)
        return rows

    header_line = lines[start].strip()
    if not header_line:
        start += 1
        header_line = lines[start].strip() if start < len(lines) else ""
    headers = header_line.split()
    if headers and headers[0].lower() in ("div", "divergence"):
        headers = headers[1:]
    bucket_by_col = [collapse_class(h) for h in headers]

    for line in lines[start + 1 :]:
        line = line.strip()
        if not line:
            continue
        fields = line.split()
        try:
            kimura_bin = int(float(fields[0]))
        except ValueError:
            continue
        if kimura_bin > landscape_max_div:
            continue
        values = fields[1:]
        for col_idx, bucket in enumerate(bucket_by_col):
            if col_idx >= len(values):
                continue
            try:
                bp = int(float(values[col_idx]))
            except ValueError:
                continue
            if bp == 0:
                continue
            rows.append((kimura_bin, bucket, bp))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-file", required=True, help="RepeatMasker .out")
    ap.add_argument("--tbl-file", required=True, help="RepeatMasker .tbl")
    ap.add_argument("--divsum-file", required=True)
    ap.add_argument("--assembly-stats", required=True)
    ap.add_argument("--arm", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--tissue", required=True)
    ap.add_argument("--landscape-max-div", type=int, default=50)
    ap.add_argument("--class-out", required=True)
    ap.add_argument("--family-out", required=True)
    ap.add_argument("--divergence-out", required=True)
    args = ap.parse_args()

    hits = list(parse_out_file(args.out_file))
    total_bp, non_n_bp = read_assembly_stats(args.assembly_stats)

    class_bp = bp_per_group(hits, collapse_class)
    family_bp = bp_per_group(hits, lambda cf: cf)

    for bucket in CANONICAL_CLASSES:
        class_bp.setdefault(bucket, 0)

    computed_total = total_masked_bp(hits)
    tbl_total = parse_tbl_masked_bp(args.tbl_file)
    if tbl_total is not None and tbl_total > 0:
        rel_diff = abs(computed_total - tbl_total) / tbl_total
        if rel_diff > 0.005:
            print(
                f"[summarize_rm] WARNING: computed masked bp ({computed_total}) differs "
                f"from .tbl bases-masked ({tbl_total}) by {rel_diff:.2%} for "
                f"{args.arm}/{args.species}",
                file=sys.stderr,
            )

    with open(args.class_out, "w") as fh:
        fh.write("arm\tspecies\ttissue\tclass\tbp\tpct_total\tpct_non_n\n")
        for cls, bp in sorted(class_bp.items()):
            pct_total = 100.0 * bp / total_bp if total_bp else 0.0
            pct_non_n = 100.0 * bp / non_n_bp if non_n_bp else 0.0
            fh.write(f"{args.arm}\t{args.species}\t{args.tissue}\t{cls}\t{bp}\t{pct_total:.4f}\t{pct_non_n:.4f}\n")

    with open(args.family_out, "w") as fh:
        fh.write("arm\tspecies\ttissue\tclass_family\tbp\tpct_total\tpct_non_n\n")
        for cf, bp in sorted(family_bp.items()):
            pct_total = 100.0 * bp / total_bp if total_bp else 0.0
            pct_non_n = 100.0 * bp / non_n_bp if non_n_bp else 0.0
            fh.write(f"{args.arm}\t{args.species}\t{args.tissue}\t{cf}\t{bp}\t{pct_total:.4f}\t{pct_non_n:.4f}\n")

    divsum_rows = parse_divsum(args.divsum_file, args.landscape_max_div)
    with open(args.divergence_out, "w") as fh:
        fh.write("arm\tspecies\tclass\tkimura_bin\tbp\n")
        for kimura_bin, bucket, bp in divsum_rows:
            fh.write(f"{args.arm}\t{args.species}\t{bucket}\t{kimura_bin}\t{bp}\n")


if __name__ == "__main__":
    main()
