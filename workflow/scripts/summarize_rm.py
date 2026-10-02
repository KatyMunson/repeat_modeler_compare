#!/usr/bin/env python3
"""Parse one (arm, species) RepeatMasker run into per-class and per-Class/Family
non-overlapping bp tables, cross-checked against the .tbl summary, plus a
long-format, overlap-resolved divergence-landscape chunk from the .align
file. Stdlib only, no pandas. Output rows already carry arm/species/tissue so the per-run
chunks this writes can be concatenated as-is by combine_summaries.py.

Non-overlapping bp: RepeatMasker's .out keeps lower-scoring hits that
overlap a better one (the lines flagged "*"), often from a different class.
Every genome base is therefore assigned to the single highest-scoring .out
hit covering it, and class / Class/Family bp are sums over those assigned
bases. Classes add up to the total masked bp, and a base is never counted
toward two classes. The total is cross-checked against the .tbl summary.

Divergence landscape: built from the .align file, not from
calcDivergenceFromAlign.pl's .divsum. The .divsum sums every alignment in
.align, and .align keeps overlapping alignments that the .out resolves
(most visibly across tandem arrays, where adjacent monomers' alignments
overlap), so its per-class bp can exceed the genome. Here each alignment
gets a Kimura divergence computed from its own transitions and
transversions (no CpG adjustment, matching the pipeline's
`calcDivergenceFromAlign.pl -noCpGMod`), and every genome base is assigned
to the single highest-scoring alignment covering it. Each base therefore
counts once, and a class's landscape sums to about its bp in the class
table (both use the same rule; checked below, with a warning past 5%).
"""

import argparse
import heapq
import math
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
    "Unknown_tandem",
}

# Unknown families that family_tandem.py found in tandem arrays are reported
# as their own class, set from --tandem-table (see relabel()).
TANDEM_CLASS = "Unknown_tandem"


# Classes left out of the divergence landscape: Kimura divergence from a
# consensus means nothing for (CA)n or AT-rich stretches. Their alignments
# still take part in overlap resolution, so the bases they cover aren't
# handed to a lower-scoring TE alignment instead.
NO_LANDSCAPE_CLASSES = {"Simple_repeat", "Low_complexity"}

_warned_unmapped = set()


def collapse_class(raw_class_family):
    class_part = raw_class_family.split("/", 1)[0]
    if class_part in CANONICAL_CLASSES:
        return class_part
    if class_part == "Other":  # e.g. Other/host_gene from reclassify_unknown.py
        return "Other"
    if raw_class_family not in _warned_unmapped:
        _warned_unmapped.add(raw_class_family)
        print(
            f"[summarize_rm] WARNING: unmapped RepeatMasker class '{raw_class_family}' -> Other",
            file=sys.stderr,
        )
    return "Other"


def relabel(family, class_family, tandem, reclass=None):
    """Class/Family for a hit of `family`. Only Unknown families change:
    array evidence first (Unknown_tandem, from family_tandem.py), then the
    reclassification table (reclassify_unknown.py: domain / Rfam / host
    protein). Returns (class_family, reclassified_from or "")."""
    if class_family != "Unknown":
        return class_family, ""
    if family in tandem:
        return TANDEM_CLASS, ""
    if reclass and family in reclass:
        return reclass[family], "Unknown"
    return class_family, ""


def parse_out_file(path, tandem=frozenset(), reclass=None):
    """Yield (query_seq, begin, end, score, (class_family, reclassified_from))
    for every hit line."""
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
            class_family = relabel(fields[9], fields[10], tandem, reclass)
            if begin > end:
                begin, end = end, begin
            yield query_seq, begin, end, int(fields[0]), class_family


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


def owned_segments(records):
    """records: iterable of (contig, begin, end, score, payload), 1-based
    inclusive. Assigns every covered base to the highest-scoring record
    covering it (ties: earliest record) and yields (bp, payload) for each
    maximal run of bases owned by one record."""
    by_contig = {}
    for idx, (contig, begin, end, score, payload) in enumerate(records):
        by_contig.setdefault(contig, []).append((begin, end, score, idx, payload))
    for recs in by_contig.values():
        # Sweep over elementary segments: at each boundary, add records that
        # start there; the top of a max-heap on score (lazily dropping
        # records that have ended) owns the segment up to the next boundary.
        recs.sort()
        boundaries = sorted({r[0] for r in recs} | {r[1] + 1 for r in recs})
        heap = []
        next_rec = 0
        for seg_start, seg_next in zip(boundaries, boundaries[1:]):
            while next_rec < len(recs) and recs[next_rec][0] <= seg_start:
                begin, end, score, idx, payload = recs[next_rec]
                heapq.heappush(heap, (-score, idx, end, payload))
                next_rec += 1
            while heap and heap[0][2] < seg_start:
                heapq.heappop(heap)
            if heap:
                yield seg_next - seg_start, heap[0][3]


def owned_bp_by_family(hits):
    """hits: (contig, begin, end, score, (class_family, reclassified_from)).
    Returns ({class_family: bp}, {class_family: bp reclassified from
    Unknown}), each base counted once, for its best hit."""
    bp = {}
    from_unknown = {}
    for seg_bp, (class_family, reclassified_from) in owned_segments(hits):
        bp[class_family] = bp.get(class_family, 0) + seg_bp
        if reclassified_from:
            from_unknown[class_family] = from_unknown.get(class_family, 0) + seg_bp
    return bp, from_unknown


def total_masked_bp(hits):
    by_contig = {}
    for query_seq, begin, end, _score, _class_family in hits:
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


def kimura_2p(transitions, transversions, length):
    """Kimura two-parameter distance, in percent. None if undefined
    (no well-characterized bases, or too diverged for the log terms)."""
    if length <= 0:
        return None
    p = transitions / length
    q = transversions / length
    a = 1 - 2 * p - q
    b = 1 - 2 * q
    if a <= 0 or b <= 0:
        return None
    return max(0.0, 100.0 * (-0.5 * math.log(a) - 0.25 * math.log(b)))


def _is_align_header(fields):
    # "score %div %del %ins query begin end (left) [C] repeat#class ..."
    if len(fields) < 9 or not fields[0].isdigit():
        return False
    try:
        float(fields[1])
        int(fields[5])
        int(fields[6])
    except ValueError:
        return False
    return fields[7].startswith("(")


def _is_seq_line(fields):
    # "  query  begin  SEQ  end" or "C repeat#class  begin  SEQ  end"
    return len(fields) >= 4 and fields[-1].isdigit() and fields[-3].isdigit()


MARKUP_CHARS = set("iv-?")


def _is_markup_line(fields):
    # The line between query and consensus: "i" transition, "v"
    # transversion, "-" gap, "?" ambiguous; blank for an identical block.
    return all(set(tok) <= MARKUP_CHARS for tok in fields)


def _acgt_count(seq):
    return sum(seq.count(c) for c in "ACGTacgt")


def _well_characterized(query_seq, subject_seq):
    """Columns where both sides are A/C/G/T. str.count only, so a
    multi-Gb .align stays fast; a column non-ACGT on both sides (N against
    a gap) is subtracted twice, which is rare enough not to matter."""
    return len(query_seq) - (len(query_seq) - _acgt_count(query_seq)) - (len(subject_seq) - _acgt_count(subject_seq))


def parse_align_file(path, tandem=frozenset(), reclass=None):
    """Yield (contig, begin, end, score, class_family, kimura_pct_or_None)
    for every alignment in a RepeatMasker .align file (-a output)."""
    try:
        fh = open(path)
    except FileNotFoundError:
        print(f"[summarize_rm] WARNING: no .align file at {path}, skipping landscape", file=sys.stderr)
        return
    header = None
    pending_query = None
    ts = tv = length = 0
    missing_class_warned = False

    def finish():
        contig, begin, end, score, class_family = header
        return contig, begin, end, score, class_family, kimura_2p(ts, tv, length)

    with fh:
        for line in fh:
            fields = line.split()
            if not fields:
                continue
            if _is_align_header(fields):
                if header is not None:
                    yield finish()
                class_family = None
                for tok in fields[8:]:
                    if "#" in tok:
                        family, class_family = tok.split("#", 1)
                        class_family = relabel(family, class_family, tandem, reclass)[0]
                        break
                if class_family is None:
                    if not missing_class_warned:
                        missing_class_warned = True
                        print(
                            f"[summarize_rm] WARNING: .align header without a repeat#class "
                            f"token, counted as Unknown: {line.strip()}",
                            file=sys.stderr,
                        )
                    class_family = "Unknown"
                begin, end = int(fields[5]), int(fields[6])
                if begin > end:
                    begin, end = end, begin
                header = (fields[4], begin, end, int(fields[0]), class_family)
                pending_query = None
                ts = tv = length = 0
                continue
            if header is None:
                continue
            if _is_seq_line(fields):
                if pending_query is None:
                    pending_query = fields[-2]
                else:
                    length += _well_characterized(pending_query, fields[-2])
                    pending_query = None
            elif pending_query is not None and _is_markup_line(fields):
                markup = "".join(fields)
                ts += markup.count("i")
                tv += markup.count("v")
    if header is not None:
        yield finish()


def resolve_landscape(alignments, landscape_max_div):
    """alignments: iterable of (contig, begin, end, score, class_family,
    kimura). Assigns every covered base to the highest-scoring alignment
    covering it (same rule as the class table), then sums bp per (class,
    Kimura bin). Returns ({(class, bin): bp}, {class: resolved_bp_all_bins})."""
    records = (
        (contig, begin, end, score, (collapse_class(class_family), kimura))
        for contig, begin, end, score, class_family, kimura in alignments
    )
    landscape = {}
    resolved_by_class = {}
    for seg_bp, (cls, kimura) in owned_segments(records):
        resolved_by_class[cls] = resolved_by_class.get(cls, 0) + seg_bp
        if cls in NO_LANDSCAPE_CLASSES or kimura is None:
            continue
        kimura_bin = int(kimura)
        if kimura_bin > landscape_max_div:
            continue
        key = (cls, kimura_bin)
        landscape[key] = landscape.get(key, 0) + seg_bp
    return landscape, resolved_by_class


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-file", required=True, help="RepeatMasker .out")
    ap.add_argument("--tbl-file", required=True, help="RepeatMasker .tbl")
    ap.add_argument("--align-file", required=True, help="RepeatMasker .align (-a)")
    ap.add_argument("--assembly-stats", required=True)
    ap.add_argument("--arm", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--tissue", required=True)
    ap.add_argument("--landscape-max-div", type=int, default=50)
    ap.add_argument("--tandem-table", help="family_tandem.tsv: report its tandem Unknown families as Unknown_tandem")
    ap.add_argument("--reclass-table", help="unknown_reclassification.tsv: relabel Unknown families it classifies")
    ap.add_argument("--class-out", required=True)
    ap.add_argument("--family-out", required=True)
    ap.add_argument("--divergence-out", required=True)
    args = ap.parse_args()

    tandem = frozenset()
    if args.tandem_table:
        from family_tandem import tandem_families
        tandem = frozenset(tandem_families(args.tandem_table))
    reclass = {}
    if args.reclass_table:
        from reclassify_unknown import reclassified
        reclass = reclassified(args.reclass_table)
    hits = list(parse_out_file(args.out_file, tandem, reclass))
    total_bp, non_n_bp = read_assembly_stats(args.assembly_stats)

    family_bp, family_bp_from_unknown = owned_bp_by_family(hits)
    class_bp = {}
    for class_family, bp in family_bp.items():
        cls = collapse_class(class_family)
        class_bp[cls] = class_bp.get(cls, 0) + bp

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
        # bp_from_unknown: the part of bp held by Unknown families that the
        # reclassification table relabelled into this Class/Family.
        fh.write("arm\tspecies\ttissue\tclass_family\tbp\tpct_total\tpct_non_n\tbp_from_unknown\n")
        for cf, bp in sorted(family_bp.items()):
            pct_total = 100.0 * bp / total_bp if total_bp else 0.0
            pct_non_n = 100.0 * bp / non_n_bp if non_n_bp else 0.0
            fh.write(f"{args.arm}\t{args.species}\t{args.tissue}\t{cf}\t{bp}\t{pct_total:.4f}\t{pct_non_n:.4f}\t"
                     f"{family_bp_from_unknown.get(cf, 0)}\n")

    del hits
    landscape, resolved_by_class = resolve_landscape(
        parse_align_file(args.align_file, tandem, reclass), args.landscape_max_div
    )
    for cls, resolved_bp in sorted(resolved_by_class.items()):
        out_bp = class_bp.get(cls, 0)
        if out_bp and abs(resolved_bp - out_bp) / out_bp > 0.05:
            print(
                f"[summarize_rm] WARNING: {args.arm}/{args.species} {cls}: .align-resolved bp "
                f"({resolved_bp}) differs from .out non-overlapping bp ({out_bp}) by "
                f"{abs(resolved_bp - out_bp) / out_bp:.1%}",
                file=sys.stderr,
            )
    with open(args.divergence_out, "w") as fh:
        fh.write("arm\tspecies\tclass\tkimura_bin\tbp\tpct_non_n\n")
        for (cls, kimura_bin), bp in sorted(landscape.items()):
            pct_non_n = 100.0 * bp / non_n_bp if non_n_bp else 0.0
            fh.write(f"{args.arm}\t{args.species}\t{cls}\t{kimura_bin}\t{bp}\t{pct_non_n:.6f}\n")


if __name__ == "__main__":
    main()
