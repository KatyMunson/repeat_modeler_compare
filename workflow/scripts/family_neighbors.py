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
    ap.add_argument("--max-gap", type=int, default=100, help="bp between a copy and its neighbour")
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
    gaps = defaultdict(list)
    out_label = {}
    for spec in args.out_file:
        sample, _, path = spec.partition("=")
        by_contig = read_out(path, fam_index, fam_names)
        print(f"[family_neighbors] {sample}: {sum(len(v) for v in by_contig.values()):,} copies "
              f"on {len(by_contig):,} contigs", file=sys.stderr)
        scan(by_contig, stats, gaps, args.max_gap, args.max_overlap)
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

    # ---- neighbour table
    with open(args.neighbors_out, "w") as out:
        out.write("family\tclass_family\tn_copies\tside\tpartner\tpartner_class\trel_strand\tn\tfrac_of_family\t"
                  "frac_of_partner\tmedian_gap\n")
        for idx in sorted(stats, key=lambda k: -stats[k].n):
            st = stats[idx]
            if st.n < args.min_copies or not is_target(cls_of(idx)):
                continue
            for side in ("5p", "3p"):
                for (p, rel), n in st.side[side].most_common(args.top):
                    out.write(f"{fam_names[idx]}\t{cls_of(idx)}\t{st.n}\t{side}\t{fam_names[p]}\t{cls_of(p)}\t"
                              f"{rel}\t{n}\t{n / st.n:.3f}\t{n / max(stats[p].n, 1):.3f}\t"
                              f"{statistics.median(gaps[(idx, side, p, rel)]):.0f}\n")

    # ---- calls
    def best_call(idx):
        st = stats[idx]
        cands = []
        for p, n in st.both.items():        # flanked by P on both sides: internal region
            cands.append((n / st.n, n, "internal_of", p, "I"))
        same5 = {p: n for (p, rel), n in st.side["5p"].items() if rel == "same"}
        same3 = {p: n for (p, rel), n in st.side["3p"].items() if rel == "same"}
        for p in set(same5) & set(same3):
            # LTR: P on its 3' side (5' LTR) in some copies, on its 5' side (3' LTR) in others;
            # from P's side, this family sits at both of P's ends
            p_st = stats[p]
            p5 = p_st.side["5p"].get((idx, "same"), 0)
            p3 = p_st.side["3p"].get((idx, "same"), 0)
            frac_p = min(p5, p3) / max(p_st.n, 1)
            n = min(same5[p], same3[p])
            if frac_p >= args.min_frac and n >= args.min_pairs and st.both.get(p, 0) < n:
                cands.append((min(frac_p, (same5[p] + same3[p]) / st.n), n, "LTR_of", p, "LTR"))
        for side, part in (("3p", "5prime_of"), ("5p", "3prime_of")):
            for (p, rel), n in st.side[side].items():
                if rel == "same":
                    cands.append((n / st.n, n, part, p, part.split("_")[0]))
        cands = [c for c in cands if c[1] >= args.min_pairs and c[0] >= args.min_frac]
        # LTR/internal patterns outrank one-sided continuation at similar support
        rank = {"LTR_of": 2, "internal_of": 2, "5prime_of": 1, "3prime_of": 1}
        return max(cands, key=lambda c: (c[0] * (1.25 if rank[c[2]] == 2 else 1), c[1]), default=None)

    calls = {}
    for idx, st in stats.items():
        if st.n >= args.min_copies and is_target(cls_of(idx)):
            calls[idx] = best_call(idx)

    def resolve(idx, c, propagated):
        frac, n, pattern, p, part = c
        pcls = cls_of(p)
        source = "label"
        if not informative(cls_of(idx), pcls) and propagated and p in propagated:
            pcls, source = propagated[p], "propagated"
        if not informative(cls_of(idx), pcls):
            return None
        conf = "high" if frac >= args.high_frac and n >= 2 * args.min_pairs else "medium"
        if source == "propagated":
            conf = "medium"
        return pcls, conf, source

    resolved = {}
    for idx, c in calls.items():
        if c:
            r = resolve(idx, c, None)
            if r:
                resolved[idx] = r
    first = {idx: r[0] for idx, r in resolved.items()}
    for idx, c in calls.items():
        if c and idx not in resolved:
            r = resolve(idx, c, first)
            if r:
                resolved[idx] = r

    n_called = Counter()
    with open(args.calls_out, "w") as out, open(args.proposals_out, "w") as prop:
        out.write("family\tclass_family\tn_copies\ttandem_family\tpolyA_3p_frac\ttandem_neighbor_frac\tpattern\t"
                  "partner\tpartner_class\tsupport_copies\tsupport_frac\tproposed_class_family\tconfidence\tnote\n")
        prop.write("group\tfamily\tclass_family\tpart\tevidence\tnote\tconfidence\n")
        for idx in sorted(calls, key=lambda k: -stats[k].n):
            st, c = stats[idx], calls[idx]
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
            out.write(f"{fam}\t{cls_of(idx)}\t{st.n}\t{fam in tandem_fam}\t{pa:.3f}\t{td:.3f}\t{pattern}\t"
                      + (f"{fam_names[c[3]]}\t{cls_of(c[3])}\t{c[1]}\t{c[0]:.3f}\t" if c else "NA\tNA\t0\t0\t")
                      + (f"{r[0]}\t{r[1]}\t" if r else "NA\tnone\t")
                      + ("; ".join(notes) or ".") + "\n")
            if r:
                n_called[r[1]] += 1
                p = fam_names[c[3]]
                group = group_of.get(p, f"with-{p}")
                note = f"{pattern} {p} ({c[1]} copies, {c[0]:.2f})" + (f"; {r[2]} class" if r[2] != "label" else "")
                prop.write(f"{group}\t{fam}\t{r[0]}\t{c[4]}\tfamily_neighbors\t{note}\t{r[1]}\n")
    print(f"[family_neighbors] {len(calls)} target families tested; proposals: "
          f"{n_called['high']} high, {n_called['medium']} medium", file=sys.stderr)


if __name__ == "__main__":
    main()
