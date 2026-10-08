#!/usr/bin/env python3
"""Are long NUMT calls really mitochondrial contigs scaffolded into the
genome? Ad hoc check, not part of the workflow.

A mito contig placed into a scaffold sits between assembly gaps: runs of
N right at, or within a few kb of, both ends of the call. A genuine NUMT
is embedded in contiguous nuclear sequence. For every call in a
numts.bed / numt_hits.bed with aligned_bp >= --min-len and identity >=
--min-identity, this reports:

  - the nearest N-run (>= --min-gap bp) on each side within the largest
    --flank, its distance from the call edge and its length;
  - N-runs per side within each --flank (default 5 and 10 kb);
  - N-runs inside the call;
  - the distance to each contig end (a call at a contig end is the same
    kind of signal: the contig may be the mito contig itself, trimmed);
  - verdict: both_sides (gap or contig end within the nearest flank on
    both sides: likely a scaffolded mito contig), one_side, none.

Expected by chance: the log gives the N-run density on the scanned
contigs and the probability that a random call of the same size would
have a gap within each flank on one side / both sides, to compare with
the observed fractions.

Usage (prepped genome: the names numts.bed uses):
  python3 workflow/scripts/numt_gap_check.py \\
      --bed results_v3/numt/Mlim/numts.bed --genome results_v3/Mlim/genome/Mlim.fa \\
      --out results_v3/numt/Mlim/numt_gap_check.tsv
Stdlib only."""

import argparse
import bisect
import math
import re
import sys

from fasta_utils import iter_fasta, seq_id


def read_calls(path, min_len, min_identity):
    calls = []
    with open(path) as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            r["start"], r["end"] = int(r["start"]), int(r["end"])
            if int(r["aligned_bp"]) >= min_len and float(r["identity"]) >= min_identity:
                calls.append(r)
    return calls


def n_runs(seq, min_gap):
    """[(start, end)] 0-based half-open runs of N/n >= min_gap."""
    return [(m.start(), m.end()) for m in re.finditer(r"[Nn]{%d,}" % min_gap, seq)]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bed", required=True, help="numts.bed or numt_hits.bed from numt_calls")
    ap.add_argument("--genome", required=True, help="the prepped genome FASTA (same contig names as the BED)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-len", type=int, default=5000, help="aligned bp a call needs to be checked")
    ap.add_argument("--min-identity", type=float, default=0.0)
    ap.add_argument("--min-gap", type=int, default=10, help="shortest N-run counted as an assembly gap")
    ap.add_argument("--flank", type=int, nargs="+", default=[5000, 10000], help="flank sizes, bp")
    args = ap.parse_args()
    flanks = sorted(args.flank)
    near, far = flanks[0], flanks[-1]

    calls = read_calls(args.bed, args.min_len, args.min_identity)
    by_contig = {}
    for c in calls:
        by_contig.setdefault(c["contig"], []).append(c)
    print(f"[gap_check] {len(calls)} calls with aligned_bp >= {args.min_len} and identity >= "
          f"{args.min_identity} on {len(by_contig)} contigs")
    if not calls:
        sys.exit("[gap_check] nothing to check")

    rows, seen = [], set()
    scanned_bp = scanned_gaps = 0
    for header, seq in iter_fasta(args.genome):
        contig = seq_id(header)
        if contig not in by_contig:
            continue
        seen.add(contig)
        L = len(seq)
        gaps = n_runs(seq, args.min_gap)
        scanned_bp += L
        scanned_gaps += len(gaps)
        g_starts = [g[0] for g in gaps]
        g_ends = [g[1] for g in gaps]
        for c in by_contig[contig]:
            s, e = c["start"], c["end"]
            # nearest gap ending at or before s (left) / starting at or after e (right)
            i = bisect.bisect_right(g_ends, s) - 1
            left = gaps[i] if i >= 0 else None
            j = bisect.bisect_left(g_starts, e)
            right = gaps[j] if j < len(gaps) else None
            left_d = s - left[1] if left else None
            right_d = right[0] - e if right else None
            inside = [g for g in gaps if g[0] < e and g[1] > s]
            row = {
                "id": c["id"], "contig": contig, "start": s, "end": e, "aligned_bp": c["aligned_bp"],
                "identity": c["identity"], "strand": c["strand"], "n_hits": c.get("n_hits", "1"),
                "contig_len": L, "dist_contig_start": s, "dist_contig_end": L - e,
                "left_gap_dist": left_d if left_d is not None and left_d <= far else "NA",
                "left_gap_len": left[1] - left[0] if left_d is not None and left_d <= far else "NA",
                "right_gap_dist": right_d if right_d is not None and right_d <= far else "NA",
                "right_gap_len": right[1] - right[0] if right_d is not None and right_d <= far else "NA",
                "gaps_inside": len(inside), "gap_bp_inside": sum(min(e, b) - max(s, a) for a, b in inside),
            }
            for f in flanks:
                row[f"gaps_left_{f}"] = sum(1 for a, b in gaps if b <= s and b > s - f)
                row[f"gaps_right_{f}"] = sum(1 for a, b in gaps if a >= e and a < e + f)
            left_hit = (left_d is not None and left_d <= near) or s <= near
            right_hit = (right_d is not None and right_d <= near) or L - e <= near
            row["verdict"] = "both_sides" if left_hit and right_hit else "one_side" if left_hit or right_hit else "none"
            rows.append(row)

    missing = set(by_contig) - seen
    if missing:
        sys.exit(f"[gap_check] contigs in the BED but not in the genome: {sorted(missing)[:5]} ... "
                 f"(use the prepped genome, whose names the BED uses)")

    cols = list(rows[0])
    rows.sort(key=lambda r: (r["contig"], r["start"]))
    with open(args.out, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join(str(r[c]) for c in cols) + "\n")

    # observed vs chance
    rate = scanned_gaps / scanned_bp if scanned_bp else 0.0
    print(f"[gap_check] scanned contigs: {scanned_bp} bp, {scanned_gaps} N-runs >= {args.min_gap} bp "
          f"({1e6 * rate:.2f} per Mb)")
    n = len(rows)
    for f in flanks:
        p1 = 1 - math.exp(-rate * f)
        obs_one = sum(1 for r in rows if r[f"gaps_left_{f}"] or r[f"gaps_right_{f}"]) / n
        obs_both = sum(1 for r in rows if r[f"gaps_left_{f}"] and r[f"gaps_right_{f}"]) / n
        print(f"[gap_check] within {f} bp: gap on either side {obs_one:.2f} (chance {1 - (1 - p1) ** 2:.3f}), "
              f"both sides {obs_both:.2f} (chance {p1 ** 2:.4f})")
    counts = {}
    for r in rows:
        counts[r["verdict"]] = counts.get(r["verdict"], 0) + 1
    print(f"[gap_check] verdicts (gap or contig end within {near} bp): "
          + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())))
    print(f"[gap_check] wrote {args.out}")


if __name__ == "__main__":
    main()
