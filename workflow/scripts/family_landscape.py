#!/usr/bin/env python3
"""Per-family divergence landscape for one (arm, sample): which families make
up a peak in summary/divergence_landscape.tsv.

Same rules as summarize_rm.py's class landscape (overlap-resolved: every base
goes to the highest-scoring .align alignment covering it; Kimura without CpG
adjustment; 1 % bins) and the same labels (pass the same --tandem-table /
--reclass-table / --curated-table the summarize rule used, see its log or the
Snakefile), so a family's bins add up to its share of the class peak.

Output (long): family, class_family, class, kimura_bin, bp, for the classes
in --classes. With --windows, also prints, per window, the families holding
the most bp in it and their share of the class bp there.

    python3 workflow/scripts/family_landscape.py \\
        --align-file OUT/shared/ind6/repeatmasker/ind6.fa.align \\
        --tandem-table OUT/shared/ind6/summary/family_tandem.tsv \\
        --curated-table OUT/summary/satellite_applied.tsv \\
        --crosscheck-table OUT/summary/satellite_applied.tsv \\
        --classes Satellite --windows 12-14 15-19 \\
        --out ind6.satellite_family_landscape.tsv
"""

import argparse

from summarize_rm import collapse_class, owned_segments, parse_align_file


def family_landscape(alignments, classes, max_div):
    """alignments: parse_align_file(..., with_family=True). Returns
    {(family, class_family, kimura_bin): bp} for the classes kept."""
    landscape = {}
    for seg_bp, ((family, class_family), kimura) in owned_segments(
        (contig, begin, end, score, (fam_cf, kimura))
        for contig, begin, end, score, fam_cf, kimura in alignments
    ):
        if kimura is None or collapse_class(class_family) not in classes:
            continue
        kimura_bin = int(kimura)
        if kimura_bin > max_div:
            continue
        key = (family, class_family, kimura_bin)
        landscape[key] = landscape.get(key, 0) + seg_bp
    return landscape


def parse_window(text):
    lo, hi = text.split("-")
    return int(lo), int(hi)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--align-file", required=True, help="RepeatMasker .align (-a)")
    ap.add_argument("--tandem-table", help="family_tandem.tsv (as summarize)")
    ap.add_argument("--reclass-table", help="unknown_reclassification.tsv (as summarize)")
    ap.add_argument("--curated-table", nargs="*", default=[], help="as summarize: first table listing a family wins")
    ap.add_argument("--crosscheck-table", nargs="*", default=[], help="as summarize")
    ap.add_argument("--classes", nargs="+", default=["Satellite"], help="classes to report (default Satellite)")
    ap.add_argument("--max-div", type=int, default=50)
    ap.add_argument("--windows", nargs="*", default=[], help="Kimura bin ranges, inclusive, e.g. 12-14 15-19")
    ap.add_argument("--top", type=int, default=10, help="families listed per window")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    tandem = frozenset()
    if args.tandem_table:
        from family_tandem import tandem_families
        tandem = frozenset(tandem_families(args.tandem_table))
    reclass = {}
    if args.reclass_table:
        from reclassify_unknown import reclassified
        reclass = reclassified(args.reclass_table)
    curated = {}
    if args.curated_table:
        from family_groups import curated_classes
        for path in args.curated_table:
            origin = "crosscheck" if path in args.crosscheck_table else "curated"
            for fam, cf in curated_classes(path).items():
                curated.setdefault(fam, (cf, origin))

    landscape = family_landscape(
        parse_align_file(args.align_file, tandem, reclass, curated, with_family=True),
        set(args.classes), args.max_div,
    )
    with open(args.out, "w") as fh:
        fh.write("family\tclass_family\tclass\tkimura_bin\tbp\n")
        for (family, cf, kimura_bin), bp in sorted(landscape.items(), key=lambda kv: (kv[0][2], -kv[1])):
            fh.write(f"{family}\t{cf}\t{collapse_class(cf)}\t{kimura_bin}\t{bp}\n")

    for window in args.windows:
        lo, hi = parse_window(window)
        by_family = {}
        for (family, cf, kimura_bin), bp in landscape.items():
            if lo <= kimura_bin <= hi:
                by_family[(family, cf)] = by_family.get((family, cf), 0) + bp
        total = sum(by_family.values())
        print(f"# Kimura {lo}-{hi} %: {total / 1e6:.2f} Mb in {', '.join(args.classes)}, {len(by_family)} families")
        for (family, cf), bp in sorted(by_family.items(), key=lambda kv: -kv[1])[: args.top]:
            print(f"{family}\t{cf}\t{bp / 1e6:.2f} Mb\t{100 * bp / total:.1f} %")


if __name__ == "__main__":
    main()
