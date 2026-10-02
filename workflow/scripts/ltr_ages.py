#!/usr/bin/env python3
"""Insertion ages of one LTR element from its intact copies: the
5'/3' LTR similarity of every LTRharvest / LTR_FINDER candidate whose
internal region is the element's. Diagnostic only, not run by the
pipeline. Stdlib only.

The two LTRs of a copy are identical when it inserts and diverge
independently afterwards, so their distance d dates the insertion:
T = d / (2 * rate). That clock is per copy and does not depend on how well
a library consensus fits, unlike the .out divergence of each hit.

  python3 workflow/scripts/ltr_ages.py \\
      --scn results_v3/Esto/ltr/rawLTR.scn \\
      --genome results_v3/Esto/genome/Esto.fa \\
      --out-file results_v3/shared/Esto/repeatmasker/Esto.fa.out \\
      --family Esto_rnd-1_family-332 --outdir diag/fam332_ages [--rate 2.2e-9]

Inputs:
  --scn     the pipeline's {species}/ltr/rawLTR.scn: 11 LTRharvest columns
            s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) s(rLTR) e(rLTR)
            l(rLTR) sim seq-nr, seq-nr = 0-based index of the contig in
            --genome (normalize_scn.py). Candidates found by both tools are
            kept once (reciprocal overlap >= 0.8, the longer one).
  --family  one or more .out families making up the element's INTERNAL
            region (e.g. the RepeatModeler family plus related pieces).
A candidate is a copy of the element when those families' .out hits cover
>= --min-internal-frac of its internal region (between the LTRs) and
>= --min-internal-bp bp.

Writes:
  copies.tsv    contig, element/LTR coordinates, LTR lengths, sim, internal
                coverage by the families, Jukes-Cantor distance,
                age (with --rate)
  summary.txt   counts, LTR length and similarity distribution, histogram in
                0.5% similarity bins, ages at the median and quartiles
"""

import argparse
import math
import os
import statistics
import sys
from collections import defaultdict

from fasta_utils import open_maybe_gz


def contig_names(genome):
    """Contig names in FASTA order (.fai if present, else the headers)."""
    fai = genome + ".fai"
    if os.path.exists(fai):
        with open(fai) as fh:
            return [line.split("\t", 1)[0] for line in fh if line.strip()]
    names = []
    with open_maybe_gz(genome) as fh:
        for line in fh:
            if line.startswith(">"):
                names.append(line[1:].split()[0])
    return names


def read_scn(path, names):
    """[(contig, s, e, ls, le, rs, re, sim)] from an 11/12-column .scn."""
    out = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            if len(f) < 11:
                continue
            try:
                s, e, _l, ls, le, _ll, rs, re_, _rl = (int(x) for x in f[:9])
                sim = float(f[9])
                idx = int(f[10])
            except ValueError:
                continue
            contig = f[11] if len(f) > 11 else names[idx]
            out.append((contig, s, e, ls, le, rs, re_, sim))
    return out


def dedupe(cands, min_recip=0.8):
    by_contig = defaultdict(list)
    for c in cands:
        by_contig[c[0]].append(c)
    kept = []
    for contig, cs in by_contig.items():
        cs.sort(key=lambda c: (c[1], -(c[2] - c[1])))
        group = []
        for c in cs:
            if group:
                g = group[-1]
                ov = min(g[2], c[2]) - max(g[1], c[1]) + 1
                if ov > 0 and ov >= min_recip * (c[2] - c[1] + 1) and ov >= min_recip * (g[2] - g[1] + 1):
                    if c[2] - c[1] > g[2] - g[1]:
                        group[-1] = c
                    continue
            group.append(c)
        kept.extend(group)
    return kept


def read_family_hits(path, families):
    """{contig: sorted [(begin, end)]} of the families' .out hits."""
    hits = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 11 or not f[0].isdigit() or f[9] not in families:
                continue
            b, e = sorted((int(f[5]), int(f[6])))
            hits[f[4]].append((b, e))
    for c in hits:
        hits[c].sort()
    return hits


def covered(intervals, s, e):
    """bp of [s, e] covered by sorted, possibly overlapping intervals."""
    import bisect
    i = bisect.bisect_left(intervals, (s - 200000, 0))
    cov, cur = 0, s
    for b, en in intervals[i:]:
        if b > e:
            break
        if en < cur:
            continue
        b = max(b, cur)
        en = min(en, e)
        if en >= b:
            cov += en - b + 1
            cur = en + 1
    return cov


def jc(p_diff):
    if p_diff >= 0.75:
        return float("inf")
    return -0.75 * math.log(1 - 4.0 * p_diff / 3.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scn", required=True)
    ap.add_argument("--genome", required=True, help="the genome the .scn indexes (contig order)")
    ap.add_argument("--out-file", required=True, help="RepeatMasker .out")
    ap.add_argument("--family", nargs="+", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-internal-frac", type=float, default=0.5)
    ap.add_argument("--min-internal-bp", type=int, default=1000)
    ap.add_argument("--rate", type=float, default=0.0, help="substitutions/site/year for ages (optional)")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    names = contig_names(args.genome)
    raw = read_scn(args.scn, names)
    cands = dedupe(raw)
    hits = read_family_hits(args.out_file, set(args.family))
    copies = []
    for contig, s, e, ls, le, rs, re_, sim in cands:
        if contig not in hits:
            continue
        int_s, int_e = le + 1, rs - 1
        if int_e <= int_s:
            continue
        cov = covered(hits[contig], int_s, int_e)
        frac = cov / (int_e - int_s + 1)
        if cov >= args.min_internal_bp and frac >= args.min_internal_frac:
            d = jc(1 - sim / 100.0)
            age = d / (2 * args.rate) if args.rate else float("nan")
            copies.append((contig, s, e, ls, le, rs, re_, sim, cov, frac, d, age))

    with open(os.path.join(args.outdir, "copies.tsv"), "w") as out:
        out.write("contig\tstart\tend\tlLTR_start\tlLTR_end\trLTR_start\trLTR_end\tlLTR_len\trLTR_len\t"
                  "ltr_sim_pct\tinternal_bp_by_family\tinternal_frac\tjc_distance\tage_years\n")
        for c in sorted(copies, key=lambda c: -c[7]):
            contig, s, e, ls, le, rs, re_, sim, cov, frac, d, age = c
            out.write(f"{contig}\t{s}\t{e}\t{ls}\t{le}\t{rs}\t{re_}\t{le - ls + 1}\t{re_ - rs + 1}\t"
                      f"{sim:.2f}\t{cov}\t{frac:.2f}\t{d:.5f}\t{age:.0f}\n")

    lines = [f"families\t{' '.join(args.family)}",
             f"scn_candidates\t{len(raw)} ({len(cands)} after merging duplicates)",
             f"intact_copies\t{len(copies)} (internal region >= {args.min_internal_frac:.0%} and "
             f">= {args.min_internal_bp} bp covered by the families)"]
    if copies:
        sims = sorted(c[7] for c in copies)
        q = statistics.quantiles(sims, n=4) if len(sims) > 1 else [sims[0]] * 3
        lens = [c[4] - c[3] + 1 for c in copies] + [c[6] - c[5] + 1 for c in copies]
        elen = [c[2] - c[1] + 1 for c in copies]
        lines += [
            f"element_len\tmedian {statistics.median(elen):.0f} (range {min(elen)}-{max(elen)})",
            f"ltr_len\tmedian {statistics.median(lens):.0f} (range {min(lens)}-{max(lens)})",
            f"ltr_sim_pct\tmedian {statistics.median(sims):.2f}; quartiles {q[0]:.2f} / {q[2]:.2f}; "
            f"min {sims[0]:.2f}; max {sims[-1]:.2f}",
            f"identical_ltrs\t{sum(1 for s in sims if s >= 99.95)} copies at 100%",
        ]
        if args.rate:
            for label, s in (("median", statistics.median(sims)), ("youngest_quartile", q[2]),
                             ("oldest_quartile", q[0])):
                lines.append(f"age_{label}\t{jc(1 - s / 100) / (2 * args.rate) / 1e6:.2f} My "
                             f"at {args.rate:g} subst/site/yr")
        hist = defaultdict(int)
        for s in sims:
            hist[math.floor(s * 2) / 2] += 1
        lines.append("sim_histogram (bin start %: copies)")
        for b in sorted(hist, reverse=True):
            lines.append(f"  {b:5.1f}\t{hist[b]}\t{'#' * min(80, hist[b])}")
    with open(os.path.join(args.outdir, "summary.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
