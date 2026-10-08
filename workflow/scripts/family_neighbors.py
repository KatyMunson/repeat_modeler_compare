#!/usr/bin/env python3
"""Place Unknown library families by what sits next to their copies.

A family that is a piece of a larger element (an LTR, the 5' end of a LINE,
the internal region of an LTR element whose LTRs are classified) has its
copies flanked, again and again, by the element's other pieces, on a fixed
side and in the same orientation. Inserted elsewhere at random, a family's
neighbours vary and their strand is 50/50. This script counts, from the
RepeatMasker .out alone, each family's neighbours on its 5' and 3' sides (in
the family's own orientation), and calls the patterns:

  LTR_of        copies have partner P (same strand) on BOTH sides, in
                different copies: 3' neighbour P = 5' LTR, 5' neighbour P =
                3' LTR. Checked from P's side too: >= --min-frac of P's
                copies carry this family on each side.
  internal_of   copies are flanked by P on both sides, same strand (the
                family is the internal region, P the LTR)
  5prime_of / 3prime_of
                copies continue into P on one side, same strand (the family
                is the 5' / 3' part of P's element: a LINE 5' UTR, a split
                consensus)
  polyA         an A-rich simple repeat right after the 3' end (non-LTR
                retrotransposon or SINE; note only)
  tandem        copies next to copies of the same family (note only)

Windows: neighbours are looked for within --max-gap bp of each copy, in
widening windows (default 100, 500, 2000). Each family takes its call from
the tightest window that gives one (pieces of young elements abut, while
old copies are split by unmasked, diverged stretches); a wider window that
finds the same pattern and partner with more support then supplies the
support (and its max_gap is reported). A wide window in a
repeat-dense genome almost always finds some neighbour, so every pattern
also needs the partner on the same strand >= --strand-ratio times as often
as on the opposite one (random neighbours are 50/50), and calls from
windows wider than --high-max-gap are capped at medium.

Targets are families labelled Unknown or <Order>/Unknown after the curated,
cross-check and reclassification labels; partners are any non-simple
family. A call needs >= --min-pairs copies and >= --min-frac of the target's
copies (and, for LTR_of, of the partner's). A partner that is itself Unknown
gives no class; one round of propagation lets a target linked to a family
called in the first round inherit that call (capped at medium).

Report-only. Proposals use the classify.curated_families format: review and
copy rows into your curated table. Stdlib only.

  python3 workflow/scripts/family_neighbors.py \\
      --out-file Esto=results_v3/shared/Esto/repeatmasker/Esto.fa.out \\
      --out-file Mlim=results_v3/shared/Mlim/repeatmasker/Mlim.fa.out \\
      --family-tandem results_v3/summary_shared_only/family_tandem.tsv \\
      --reclass results_v3/classify/unknown_reclassification.tsv \\
      --curated curation/hagfish_curated.tsv results_v3/summary/satellite_applied.tsv \\
      --neighbors-out neighbors.tsv --calls-out calls.tsv --proposals-out proposals.tsv
"""

import argparse
import bisect
import statistics
import sys
from collections import Counter, defaultdict

from family_groups import curated_classes, read_groups

SIMPLE_CLASSES = ("Simple_repeat", "Low_complexity")
COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def bare(name):
    return name.split("#", 1)[0]


def a_rich(unit, min_a=0.75):
    unit = unit.upper()
    return bool(unit) and unit.count("A") / len(unit) >= min_a


def simple_unit(name):
    """'(CA)n' -> 'CA'; 'A-rich' -> 'A'; else ''."""
    if name.startswith("(") and ")n" in name:
        return name[1:name.index(")")]
    if name.endswith("-rich"):
        return name.split("-", 1)[0]
    return ""


def read_out(path, fam_index, fam_names):
    """{contig: [(start, end, fam_idx, strand, is_simple, unit_genome)]}, one
    entry per RepeatMasker copy (hits sharing contig, ID and family merged)."""
    copies = {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 15 or not f[0].isdigit():
                continue
            contig, b, e = f[4], int(f[5]), int(f[6])
            strand = "+" if f[8] == "+" else "-"
            fam, cls = bare(f[9]), f[10]
            key = (contig, f[14], fam)
            c = copies.get(key)
            if c:
                c[0] = min(c[0], b)
                c[1] = max(c[1], e)
                continue
            if fam not in fam_index:
                fam_index[fam] = len(fam_names)
                fam_names.append(fam)
            is_simple = cls.split("/", 1)[0] in SIMPLE_CLASSES
            unit = simple_unit(fam) if is_simple else ""
            if unit and strand == "-":
                unit = unit.translate(COMP)[::-1]
            copies[key] = [b, e, fam_index[fam], strand, is_simple, unit]
    by_contig = defaultdict(list)
    for (contig, _id, _fam), c in copies.items():
        by_contig[contig].append(tuple(c))
    for contig in by_contig:
        by_contig[contig].sort()
    return by_contig


class Stats:
    __slots__ = ("n", "side", "both", "polya", "tandem")

    def __init__(self):
        self.n = 0
        self.side = {"5p": Counter(), "3p": Counter()}   # (partner_idx, rel_strand) -> copies
        self.both = Counter()                            # partner_idx -> copies flanked same-strand both sides
        self.polya = 0
        self.tandem = 0


def scan(by_contig, stats, gaps, max_gap, max_overlap):
    for contig, cs in by_contig.items():
        starts = [c[0] for c in cs]
        by_end = sorted(range(len(cs)), key=lambda k: cs[k][1])
        ends = [cs[k][1] for k in by_end]
        for i, (s, e, fam, strand, is_simple, _u) in enumerate(cs):
            if is_simple:
                continue
            # genome-downstream: starts in (e - max_overlap, e + max_gap], ends past e
            down = down_simple = None
            j = bisect.bisect_right(starts, s)
            while j < len(cs) and cs[j][0] <= e + max_gap:
                c = cs[j]
                if j != i and c[1] > e and c[0] > s and c[0] >= e - max_overlap:
                    if c[4]:
                        if down_simple is None:
                            down_simple = c
                    elif down is None or c[0] < down[0]:
                        down = c
                j += 1
            # genome-upstream: ends in [s - max_gap, s + max_overlap), starts before s
            up = up_simple = None
            lo = bisect.bisect_left(ends, s - max_gap)
            hi = bisect.bisect_left(ends, s + max_overlap)
            for k in range(lo, hi):
                c = cs[by_end[k]]
                if by_end[k] != i and c[0] < s and c[1] < e:
                    if c[4]:
                        if up_simple is None or c[1] > up_simple[1]:
                            up_simple = c
                    elif up is None or c[1] > up[1]:
                        up = c
            st = stats[fam]
            st.n += 1
            # family orientation: + strand -> 5' is genome-upstream
            five, three = (up, down) if strand == "+" else (down, up)
            three_simple = down_simple if strand == "+" else up_simple
            if three_simple is not None:
                unit = three_simple[5]
                if strand == "-":
                    unit = unit.translate(COMP)[::-1]
                if a_rich(unit):
                    st.polya += 1
            sides = {}
            for side, nb in (("5p", five), ("3p", three)):
                if nb is None:
                    continue
                rel = "same" if nb[3] == strand else "opp"
                if nb[2] == fam:
                    st.tandem += 1
                    continue
                st.side[side][(nb[2], rel)] += 1
                sides[side] = (nb[2], rel)
                gap = (nb[0] - e) if nb[0] > s else (s - nb[1])
                gaps[(fam, side, nb[2], rel)].append(gap)
            if "5p" in sides and sides.get("5p") == sides.get("3p") and sides["5p"][1] == "same":
                st.both[sides["5p"][0]] += 1


def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            yield dict(zip(header, line.rstrip("\n").split("\t")))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-file", nargs="+", required=True, help="sample=path to RepeatMasker .out")
    ap.add_argument("--family-tandem", default="", help="combined family_tandem.tsv (labels, tandem flag)")
    ap.add_argument("--reclass", default="", help="unknown_reclassification.tsv")
    ap.add_argument("--curated", nargs="*", default=[], help="curated / applied tables (first wins)")
    ap.add_argument("--max-gap", type=int, nargs="+", default=[100, 500, 2000],
                    help="windows (bp between a copy and its neighbour), tried tightest first: each family "
                         "takes its call from the first window that gives one")
    ap.add_argument("--high-max-gap", type=int, default=500,
                    help="calls from wider windows are capped at medium")
    ap.add_argument("--strand-ratio", type=float, default=3.0,
                    help="same-strand partner copies >= this x opposite-strand ones (random neighbours are 50/50)")
    ap.add_argument("--max-overlap", type=int, default=50, help="bp a neighbour may overlap the copy")
    ap.add_argument("--min-copies", type=int, default=50, help="target copies (all samples) to be tested")
    ap.add_argument("--min-pairs", type=int, default=20, help="copies showing the pattern")
    ap.add_argument("--min-frac", type=float, default=0.2, help="medium: share of copies showing the pattern")
    ap.add_argument("--high-frac", type=float, default=0.4, help="high: share of copies showing the pattern")
    ap.add_argument("--top", type=int, default=5, help="partners per side in --neighbors-out")
    ap.add_argument("--neighbors-out", required=True)
    ap.add_argument("--calls-out", required=True)
    ap.add_argument("--proposals-out", required=True)
    args = ap.parse_args()

    # labels: curated (first table wins) > reclassification > family_tandem
    label, group_of, tandem_fam = {}, {}, set()
    if args.family_tandem:
        for r in read_tsv(args.family_tandem):
            label.setdefault(bare(r["family"]), r["class_family"])
            if r.get("tandem_family") == "True":
                tandem_fam.add(bare(r["family"]))
    if args.reclass:
        for r in read_tsv(args.reclass):
            if r.get("new_class") and r["new_class"] != r.get("old_class"):
                label[bare(r["family"])] = r["new_class"]
    cur = {}
    for path in args.curated:
        for fam, cf in curated_classes(path).items():
            cur.setdefault(fam, cf)
        for g, fam, _row in read_groups(path):
            group_of.setdefault(fam, g)
    label.update(cur)

    fam_index, fam_names = {}, []
    stats = defaultdict(Stats)
    tiers = sorted(set(args.max_gap))
    stats = {t: defaultdict(Stats) for t in tiers}
    gaps = {t: defaultdict(list) for t in tiers}
    out_label = {}
    for spec in args.out_file:
        sample, _, path = spec.partition("=")
        by_contig = read_out(path, fam_index, fam_names)
        print(f"[family_neighbors] {sample}: {sum(len(v) for v in by_contig.values()):,} copies "
              f"on {len(by_contig):,} contigs", file=sys.stderr)
        for t in tiers:
            scan(by_contig, stats[t], gaps[t], t, args.max_overlap)
        del by_contig
        # .out labels for families missing from the tables
        with open(path) as fh:
            for line in fh:
                f = line.split()
                if len(f) >= 15 and f[0].isdigit():
                    out_label.setdefault(bare(f[9]), f[10])
    for fam, cls in out_label.items():
        label.setdefault(fam, cls)

    def cls_of(idx):
        return label.get(fam_names[idx], "Unknown")

    def is_target(cls):
        return cls == "Unknown" or cls.endswith("/Unknown")

    def informative(target_cls, partner_cls):
        """partner label usable for the target: classified (an <Order>/Unknown
        partner only lifts a plain Unknown to that order), and compatible with
        the target's order when the target is <Order>/Unknown"""
        if partner_cls == "Unknown" or partner_cls.split("/", 1)[0] in SIMPLE_CLASSES + ("Satellite",):
            return False
        if is_target(partner_cls):           # <Order>/Unknown: gives a plain Unknown its order
            return target_cls == "Unknown"
        if target_cls != "Unknown":
            return partner_cls.split("/", 1)[0] == target_cls.split("/", 1)[0]
        return True

    base = stats[tiers[0]]      # copy counts and the poly(A) / tandem notes: tightest window
    targets = [idx for idx, st in base.items() if st.n >= args.min_copies and is_target(cls_of(idx))]
    targets.sort(key=lambda k: -base[k].n)

    # ---- neighbour table (every window)
    with open(args.neighbors_out, "w") as out:
        out.write("family\tclass_family\tn_copies\tmax_gap\tside\tpartner\tpartner_class\trel_strand\tn\t"
                  "frac_of_family\tfrac_of_partner\tmedian_gap\n")
        for idx in targets:
            for t in tiers:
                st = stats[t][idx]
                for side in ("5p", "3p"):
                    for (p, rel), n in st.side[side].most_common(args.top):
                        out.write(f"{fam_names[idx]}\t{cls_of(idx)}\t{st.n}\t{t}\t{side}\t{fam_names[p]}\t"
                                  f"{cls_of(p)}\t{rel}\t{n}\t{n / st.n:.3f}\t"
                                  f"{n / max(stats[t][p].n, 1):.3f}\t"
                                  f"{statistics.median(gaps[t][(idx, side, p, rel)]):.0f}\n")

    # ---- calls
    def biased(same, opp):
        """same-strand partner well above the 50/50 of random neighbours"""
        return same >= args.strand_ratio * max(opp, 0) and same > opp

    def best_call(idx, S):
        st = S[idx]
        cands = []
        for p, n in st.both.items():        # flanked by P on both sides: internal region
            cands.append((n / st.n, n, "internal_of", p, "I"))
        same5 = {p: n for (p, rel), n in st.side["5p"].items() if rel == "same"}
        same3 = {p: n for (p, rel), n in st.side["3p"].items() if rel == "same"}
        for p in set(same5) & set(same3):
            # LTR: P on its 3' side (5' LTR) in some copies, on its 5' side (3' LTR) in others;
            # from P's side, this family sits at both of P's ends
            if not (biased(same5[p], st.side["5p"].get((p, "opp"), 0))
                    and biased(same3[p], st.side["3p"].get((p, "opp"), 0))):
                continue
            p_st = S[p]
            p5 = p_st.side["5p"].get((idx, "same"), 0)
            p3 = p_st.side["3p"].get((idx, "same"), 0)
            frac_p = min(p5, p3) / max(p_st.n, 1)
            n = min(same5[p], same3[p])
            if frac_p >= args.min_frac and n >= args.min_pairs and st.both.get(p, 0) < n:
                cands.append((min(frac_p, (same5[p] + same3[p]) / st.n), n, "LTR_of", p, "LTR"))
        for side, part in (("3p", "5prime_of"), ("5p", "3prime_of")):
            for (p, rel), n in st.side[side].items():
                if rel == "same" and biased(n, st.side[side].get((p, "opp"), 0)):
                    cands.append((n / st.n, n, part, p, part.split("_")[0]))
        cands = [c for c in cands if c[1] >= args.min_pairs and c[0] >= args.min_frac]
        # LTR/internal patterns outrank one-sided continuation at similar support
        rank = {"LTR_of": 2, "internal_of": 2, "5prime_of": 1, "3prime_of": 1}
        return max(cands, key=lambda c: (c[0] * (1.25 if rank[c[2]] == 2 else 1), c[1]), default=None)

    def resolve(idx, c, propagated, t):
        frac, n, pattern, p, part = c
        pcls = cls_of(p)
        source = "label"
        if not informative(cls_of(idx), pcls) and propagated and p in propagated:
            pcls, source = propagated[p], "propagated"
        if not informative(cls_of(idx), pcls):
            return None
        conf = "high" if frac >= args.high_frac and n >= 2 * args.min_pairs else "medium"
        if source == "propagated" or t > args.high_max_gap:
            conf = "medium"
        return pcls, conf, source

    # Per family, the tightest window giving a usable call; failing that, the
    # tightest giving any call (a link to another Unknown, for propagation).
    per_tier = {t: {idx: best_call(idx, stats[t]) for idx in targets} for t in tiers}
    calls, resolved = {}, {}
    for idx in targets:
        calls[idx] = None
        for t in tiers:
            c = per_tier[t][idx]
            if c:
                r = resolve(idx, c, None, t)
                if r:
                    calls[idx], resolved[idx] = (t, c), r
                    break
                if calls[idx] is None:
                    calls[idx] = (t, c)
    first = {idx: r[0] for idx, r in resolved.items()}
    for idx in targets:
        if idx in resolved:
            continue
        for t in tiers:
            c = per_tier[t][idx]
            r = resolve(idx, c, first, t) if c else None
            if r:
                calls[idx], resolved[idx] = (t, c), r
                break

    # The same pattern and partner in a wider window with more support (old
    # copies split by longer unmasked stretches): report that support.
    for idx, tc in calls.items():
        if not tc or not tc[1]:
            continue
        t, c = tc
        for t2 in tiers:
            c2 = per_tier[t2][idx] if t2 > t else None
            if c2 and c2[2:4] == c[2:4] and c2[0] > c[0]:
                r2 = resolve(idx, c2, first, t2) if idx in resolved else None
                calls[idx] = (t2, c2)
                if r2:
                    resolved[idx] = r2
                t, c = t2, c2

    n_called = Counter()
    by_tier = Counter()
    with open(args.calls_out, "w") as out, open(args.proposals_out, "w") as prop:
        out.write("family\tclass_family\tn_copies\ttandem_family\tpolyA_3p_frac\ttandem_neighbor_frac\tmax_gap\t"
                  "pattern\tpartner\tpartner_class\tsupport_copies\tsupport_frac\tproposed_class_family\t"
                  "confidence\tnote\n")
        prop.write("group\tfamily\tclass_family\tpart\tevidence\tnote\tconfidence\n")
        for idx in targets:
            st = base[idx]
            t, c = calls[idx] if calls[idx] else ("NA", None)
            fam = fam_names[idx]
            notes = []
            pa, td = st.polya / st.n, st.tandem / st.n
            if pa >= args.min_frac:
                notes.append("3' poly(A): non-LTR retrotransposon or SINE?")
            if td >= args.min_frac:
                notes.append("copies next to each other (tandem)")
            r = resolved.get(idx)
            if c and not r:
                notes.append(f"linked to {fam_names[c[3]]} ({cls_of(c[3])}): classify that first")
            pattern = c[2] if c else "none"
            out.write(f"{fam}\t{cls_of(idx)}\t{st.n}\t{fam in tandem_fam}\t{pa:.3f}\t{td:.3f}\t{t}\t{pattern}\t"
                      + (f"{fam_names[c[3]]}\t{cls_of(c[3])}\t{c[1]}\t{c[0]:.3f}\t" if c else "NA\tNA\t0\t0\t")
                      + (f"{r[0]}\t{r[1]}\t" if r else "NA\tnone\t")
                      + ("; ".join(notes) or ".") + "\n")
            if r:
                n_called[r[1]] += 1
                by_tier[t] += 1
                p = fam_names[c[3]]
                group = group_of.get(p, f"with-{p}")
                note = (f"{pattern} {p} ({c[1]} copies, {c[0]:.2f}, gap <= {t} bp)"
                        + (f"; {r[2]} class" if r[2] != "label" else ""))
                prop.write(f"{group}\t{fam}\t{r[0]}\t{c[4]}\tfamily_neighbors\t{note}\t{r[1]}\n")
    print(f"[family_neighbors] {len(targets)} target families tested; proposals: "
          f"{n_called['high']} high, {n_called['medium']} medium; by window: "
          + ", ".join(f"<= {t} bp: {by_tier[t]}" for t in tiers), file=sys.stderr)


if __name__ == "__main__":
    main()
