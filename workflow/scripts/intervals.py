#!/usr/bin/env python3
"""Interval helpers shared by the region-composition scripts
(ltr_skipped_composition.py, trc_crosscheck.py). Regions are 0-based
half-open (BED-like); RepeatMasker hits are 1-based inclusive, as in
summarize_rm.owned_segments. Stdlib only."""

import bisect


def make_disjoint(intervals):
    """[(start, end, tag)] 0-based half-open -> sorted, non-overlapping list.
    Where two intervals overlap, the earlier-starting one keeps the shared
    bases and the later one is trimmed (dropped if nothing is left)."""
    out = []
    for start, end, tag in sorted(intervals, key=lambda x: (x[0], x[1])):
        if out and start < out[-1][1]:
            start = out[-1][1]
        if start < end:
            out.append((start, end, tag))
    return out


def clip_indexed(hits, by_contig):
    """Cut hits down to the parts inside by_contig's regions.

    hits: (contig, begin, end, score, payload), 1-based inclusive.
    by_contig: {contig: sorted, disjoint [(start, end, ...)]}, 0-based
    half-open; anything after `end` is ignored.
    Yields (contig, begin, end, score, payload, region_index), still 1-based
    inclusive, one piece per (hit, region) overlap."""
    starts = {c: [iv[0] for iv in ivs] for c, ivs in by_contig.items()}
    for contig, begin, end, score, payload in hits:
        intervals = by_contig.get(contig)
        if not intervals:
            continue
        b0, e0 = begin - 1, end  # half-open
        i = max(bisect.bisect_right(starts[contig], b0) - 1, 0)
        while i < len(intervals) and intervals[i][0] < e0:
            s, e = intervals[i][0], intervals[i][1]
            lo, hi = max(b0, s), min(e0, e)
            if lo < hi:
                yield contig, lo + 1, hi, score, payload, i
            i += 1


def clip(hits, by_contig):
    """clip_indexed without the region index."""
    for contig, begin, end, score, payload, _i in clip_indexed(hits, by_contig):
        yield contig, begin, end, score, payload
