#!/usr/bin/env python3
"""Cross-check one sample's shared-arm RepeatMasker annotation against its
TideCluster arrays (tidecluster_regions.py --regions).

Every base is given to its single highest-scoring .out hit, the rule
class_composition.tsv uses, so a base counts once.

--trc-out, one row per TRC plus one `ALL` row:
  n_arrays, array_bp            TideCluster arrays of the TRC
  masked_bp, unmasked_bp        array bp RepeatMasker covers / misses; the
                                ALL row's unmasked_bp is satellite sequence
                                the library misses altogether
  tandem_family_bp              part held by family_tandem.py tandem families
  top_families                  family:bp:class for the families holding the
                                most array bp (up to --top)
  (ALL only) tidehunter_outside_trc_bp: TideHunter arrays that TideCluster
                                did not cluster into any TRC
--family-out, one row per library family with bp inside TRCs or TideHunter
arrays:
  owned_bp                      genome-wide (family_tandem.tsv)
  in_trc_bp, frac_in_trc        how much of the family sits in TRC arrays
  trc_bp                        TRC:bp for each TRC holding >= 1% of in_trc_bp
  tidehunter_bp                 owned bp inside TideHunter arrays
  th_monomer_median, th_monomer_support
                                bp-weighted median TideHunter monomer length
                                over those arrays, and the share of that bp
                                whose monomer is within 5% of it
Stdlib only."""

import argparse
import sys

from intervals import clip_indexed, make_disjoint
from ltr_skipped_composition import read_out_hits
from summarize_rm import owned_segments


def read_regions(path):
    trc, th = {}, {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            c, s, e, src = f[ix["contig"]], int(f[ix["start"]]), int(f[ix["end"]]), f[ix["source"]]
            if src == "trc":
                trc.setdefault(c, []).append((s, e, f[ix["id"]]))
            else:
                ml = f[ix["monomer_len"]]
                th.setdefault(c, []).append((s, e, float(ml) if ml not in ("NA", "") else None))
    return ({c: make_disjoint(v) for c, v in trc.items()},
            {c: make_disjoint(v) for c, v in th.items()})


def read_family_table(path):
    fams = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            fams[f[ix["family"]]] = {
                "class_family": f[ix["class_family"]],
                "tandem_family": f[ix["tandem_family"]],
                "owned_bp": int(f[ix["owned_bp"]]),
            }
    return fams


def owned_by_region(hits, regions):
    """{(family, contig, region_index): owned bp}, plus class per family."""
    pieces = ((c, b, e, sc, (p, c, i)) for c, b, e, sc, p, i in clip_indexed(hits, regions))
    out, cls_of, tandem = {}, {}, {}
    for bp, ((cls, family, is_tandem), contig, i) in owned_segments(pieces):
        key = (family, contig, i)
        out[key] = out.get(key, 0) + bp
        cls_of[family] = cls
        tandem[family] = is_tandem
    return out, cls_of, tandem


def subtract_bp(a, b):
    """bp of disjoint sorted intervals `a` not covered by disjoint sorted `b`."""
    total, j = 0, 0
    for s, e, *_ in a:
        cur = s
        while j < len(b) and b[j][1] <= cur:
            j += 1
        k = j
        while k < len(b) and b[k][0] < e:
            if b[k][0] > cur:
                total += b[k][0] - cur
            cur = max(cur, b[k][1])
            k += 1
        if cur < e:
            total += e - cur
    return total


def weighted_median(pairs):
    pairs = sorted(pairs)
    total = sum(w for _, w in pairs)
    acc = 0
    for v, w in pairs:
        acc += w
        if acc >= total / 2:
            return v
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-file", required=True, help="shared-arm RepeatMasker .out")
    ap.add_argument("--tandem-table", required=True, help="shared-arm family_tandem.tsv of this sample")
    ap.add_argument("--regions", required=True, help="tidecluster_regions.py --regions")
    ap.add_argument("--species", required=True)
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("--trc-out", required=True)
    ap.add_argument("--family-out", required=True)
    args = ap.parse_args()

    fams = read_family_table(args.tandem_table)
    tandem = frozenset(f for f, v in fams.items() if v["tandem_family"] == "True")
    hits = read_out_hits(args.out_file, tandem)
    trc_regions, th_regions = read_regions(args.regions)
    print(f"[trc_crosscheck] {len(hits)} hits; {sum(map(len, trc_regions.values()))} TRC arrays, "
          f"{sum(map(len, th_regions.values()))} TideHunter arrays (de-overlapped)", file=sys.stderr)

    # --- per TRC ---
    in_trc, cls_of, is_tandem = owned_by_region(hits, trc_regions)
    trc_of = {(c, i): r[2] for c, rs in trc_regions.items() for i, r in enumerate(rs)}
    per_trc = {}
    for (family, c, i), bp in in_trc.items():
        d = per_trc.setdefault(trc_of[(c, i)], {})
        d[family] = d.get(family, 0) + bp
    arrays, array_bp = {}, {}
    for c, rs in trc_regions.items():
        for s, e, trc in rs:
            arrays[trc] = arrays.get(trc, 0) + 1
            array_bp[trc] = array_bp.get(trc, 0) + (e - s)

    th_outside = sum(subtract_bp(th_regions[c], trc_regions.get(c, [])) for c in th_regions)

    with open(args.trc_out, "w") as out:
        out.write("species\ttrc\tn_arrays\tarray_bp\tmasked_bp\tunmasked_bp\tpct_masked\t"
                  "tandem_family_bp\ttop_families\ttidehunter_outside_trc_bp\n")
        tot = {"n": 0, "bp": 0, "m": 0, "t": 0}
        for trc in sorted(array_bp, key=lambda t: int(t.split("_")[1])):
            fam_bp = per_trc.get(trc, {})
            masked = sum(fam_bp.values())
            tbp = sum(bp for f, bp in fam_bp.items() if is_tandem.get(f))
            top = sorted(fam_bp.items(), key=lambda kv: -kv[1])[: args.top]
            top_s = ";".join(f"{f}:{bp}:{cls_of[f]}" for f, bp in top) or "."
            out.write(f"{args.species}\t{trc}\t{arrays[trc]}\t{array_bp[trc]}\t{masked}\t"
                      f"{array_bp[trc] - masked}\t{100.0 * masked / array_bp[trc]:.2f}\t{tbp}\t{top_s}\tNA\n")
            tot["n"] += arrays[trc]; tot["bp"] += array_bp[trc]; tot["m"] += masked; tot["t"] += tbp
        pct = 100.0 * tot["m"] / tot["bp"] if tot["bp"] else 0.0
        out.write(f"{args.species}\tALL\t{tot['n']}\t{tot['bp']}\t{tot['m']}\t{tot['bp'] - tot['m']}\t"
                  f"{pct:.2f}\t{tot['t']}\t.\t{th_outside}\n")

    # --- per family ---
    fam_trc = {}
    for (family, c, i), bp in in_trc.items():
        d = fam_trc.setdefault(family, {})
        trc = trc_of[(c, i)]
        d[trc] = d.get(trc, 0) + bp
    in_th, cls_th, _ = owned_by_region(hits, th_regions)
    mono_of = {(c, i): r[2] for c, rs in th_regions.items() for i, r in enumerate(rs)}
    fam_th = {}
    for (family, c, i), bp in in_th.items():
        fam_th.setdefault(family, []).append((mono_of[(c, i)], bp))
        cls_of.setdefault(family, cls_th[family])

    with open(args.family_out, "w") as out:
        out.write("species\tfamily\tclass_family\tclass\ttandem_family\towned_bp\tin_trc_bp\tfrac_in_trc\t"
                  "n_trcs\ttop_trc\ttop_trc_bp\ttrc_bp\ttidehunter_bp\tth_monomer_median\tth_monomer_support\n")
        for family in sorted(set(fam_trc) | set(fam_th)):
            meta = fams.get(family, {"class_family": "NA", "tandem_family": "NA", "owned_bp": 0})
            d = fam_trc.get(family, {})
            in_bp = sum(d.values())
            owned = meta["owned_bp"]
            frac = in_bp / owned if owned else 0.0
            ordered = sorted(d.items(), key=lambda kv: -kv[1])
            keep = [f"{t}:{bp}" for t, bp in ordered if bp >= 0.01 * in_bp]
            th = [(m, bp) for m, bp in fam_th.get(family, []) if m is not None]
            th_bp = sum(bp for _, bp in fam_th.get(family, []))
            med = weighted_median(th) if th else None
            sup = (sum(bp for m, bp in th if abs(m - med) <= 0.05 * med) / sum(bp for _, bp in th)
                   if med else None)
            out.write(f"{args.species}\t{family}\t{meta['class_family']}\t{cls_of.get(family, 'NA')}\t"
                      f"{meta['tandem_family']}\t{owned}\t{in_bp}\t{frac:.3f}\t{len(d)}\t"
                      f"{ordered[0][0] if ordered else 'NA'}\t{ordered[0][1] if ordered else 0}\t"
                      f"{';'.join(keep) or '.'}\t{th_bp}\t"
                      f"{'NA' if med is None else f'{med:g}'}\t{'NA' if sup is None else f'{sup:.2f}'}\n")


if __name__ == "__main__":
    main()
