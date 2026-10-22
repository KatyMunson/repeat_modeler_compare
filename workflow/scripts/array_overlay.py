#!/usr/bin/env python3
"""How does the RepeatMasker pipeline label the bases inside an external set
of tandem arrays (e.g. the satellite pipeline's stage 03 BED)? One-off
diagnostic, not run by the workflow. Stdlib only.

Each array base gets the label summarize_rm.py would give it: the best-
scoring .out hit covering it (each base counted once), relabelled exactly
as in family_composition.tsv (curated / applied cross-check tables first,
then Unknown_tandem, then the reclassification table). Bases no hit covers
are `unmasked`.

  python3 workflow/scripts/array_overlay.py \\
      --bed stage03_tandem_arrays.bed --sample Esto \\
      --out-file results_v3/shared/Esto/repeatmasker/Esto.fa.out \\
      --tandem-table results_v3/summary_shared_only/family_tandem.tsv \\
      --reclass-table results_v3/classify/unknown_reclassification.tsv \\
      --curated-table curation/hagfish_curated.tsv results_v3/summary/satellite_applied.tsv \\
      --crosscheck-table results_v3/summary/satellite_applied.tsv \\
      --name-map results_v3/Esto/genome/Esto.name_map.tsv \\
      --outdir array_overlay/Esto

--bed: contig, start, end (0-based half-open), optional 4th column (motif or
array name, used for the per-motif table). Overlapping arrays are merged so
no base counts twice. --name-map (prep_genome's original -> sanitized names)
translates BED contig names to the .out's; BED names already in the .out are
kept as they are.

Writes to --outdir:
  by_class.tsv         array bp per collapsed class (Satellite, Unknown, LTR, ...)
  by_class_family.tsv  array bp per Class/Family label
  by_family.tsv        array bp per library family: its label, its total
                       owned bp genome-wide and the share of that inside arrays
                       (a family mostly inside arrays but not labelled
                       Satellite is a relabelling candidate)
  by_motif.tsv         with a 4th BED column: per motif, array bp and the
                       label holding most of it
  summary.txt          headline numbers
"""

import argparse
import os
import sys
from collections import Counter, defaultdict

from intervals import clip_indexed, make_disjoint
from summarize_rm import collapse_class, owned_segments, relabel


def read_tandem(path, sample):
    """Families flagged tandem_family in this sample (combined or per-sample table)."""
    names = set()
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_fam, i_flag = header.index("family"), header.index("tandem_family")
        i_s = header.index("sample") if "sample" in header else None
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[i_flag] == "True" and (i_s is None or f[i_s] == sample):
                names.add(f[i_fam])
    return frozenset(names)


def read_bed(path, name_map, known):
    regions = defaultdict(list)
    unknown_contigs = Counter()
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            contig = f[0]
            if contig not in known and contig in name_map:
                contig = name_map[contig]
            if contig not in known:
                unknown_contigs[f[0]] += 1
            regions[contig].append((int(f[1]), int(f[2]), f[3] if len(f) > 3 else "."))
    return regions, unknown_contigs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bed", required=True, help="tandem arrays, 0-based half-open")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out-file", required=True, help="RepeatMasker .out (shared arm)")
    ap.add_argument("--tandem-table", default="", help="family_tandem.tsv (combined or per-sample)")
    ap.add_argument("--reclass-table", default="")
    ap.add_argument("--curated-table", nargs="*", default=[], help="curated tables, first wins")
    ap.add_argument("--crosscheck-table", nargs="*", default=[],
                    help="which --curated-table files are the applied satellite cross-check")
    ap.add_argument("--name-map", default="", help="prep_genome name_map.tsv (original -> sanitized)")
    ap.add_argument("--top", type=int, default=50, help="rows in by_family.tsv")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    tandem = read_tandem(args.tandem_table, args.sample) if args.tandem_table else frozenset()
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

    hits = []
    with open(args.out_file) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 11 or not f[0].isdigit():
                continue
            try:
                b, e = int(f[5]), int(f[6])
            except ValueError:
                continue
            fam = f[9].split("#", 1)[0]
            cf, origin = relabel(fam, f[10], tandem, reclass, curated)
            hits.append((f[4], min(b, e), max(b, e), int(f[0]), (cf, origin, fam)))
    contigs = {h[0] for h in hits}
    print(f"[array_overlay] {args.sample}: {len(hits):,} hits on {len(contigs):,} contigs", file=sys.stderr)

    name_map = {}
    if args.name_map:
        with open(args.name_map) as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2:
                    name_map[f[0]] = f[1]
    raw, unknown_contigs = read_bed(args.bed, name_map, contigs)
    if unknown_contigs:
        n = sum(unknown_contigs.values())
        print(f"[array_overlay] WARNING {n} BED intervals on {len(unknown_contigs)} contigs with no .out hits "
              f"(name mismatch, or arrays RepeatMasker masked nothing on): "
              f"{', '.join(list(unknown_contigs)[:5])}", file=sys.stderr)
    # merge overlapping arrays so each base counts once (keeps the first name)
    by_contig, array_bp = {}, 0
    for contig, ivs in raw.items():
        merged = []
        for s, e, name in sorted(ivs):
            if merged and s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e), merged[-1][2])
            else:
                merged.append((s, e, name))
        by_contig[contig] = make_disjoint(merged)
        array_bp += sum(e - s for s, e, _n in by_contig[contig])
    names = {c: [iv[2] for iv in ivs] for c, ivs in by_contig.items()}

    # each array base: its best .out hit (owned_segments over the clipped pieces)
    clipped = [(c, b, e, sc, (pl, i)) for c, b, e, sc, pl, i in clip_indexed(hits, by_contig)]
    cf_bp, cls_bp, fam_bp, origin_bp = Counter(), Counter(), Counter(), Counter()
    fam_label = {}
    masked = 0
    for seg_bp, ((cf, origin, fam), i_region) in owned_segments(clipped):
        masked += seg_bp
        cf_bp[cf] += seg_bp
        cls_bp[collapse_class(cf)] += seg_bp
        fam_bp[fam] += seg_bp
        fam_label[fam] = (cf, origin)
        origin_bp[(collapse_class(cf), origin or "library")] += seg_bp
    unmasked = array_bp - masked
    cls_bp["unmasked"] += unmasked
    cf_bp["unmasked"] += unmasked

    # genome-wide owned bp per family, for the share inside arrays
    total_fam = Counter()
    for seg_bp, (_cf, _o, fam) in owned_segments(hits):
        total_fam[fam] += seg_bp

    def pct(x):
        return f"{100.0 * x / array_bp:.2f}" if array_bp else "NA"

    with open(os.path.join(args.outdir, "by_class.tsv"), "w") as out:
        out.write("sample\tclass\tbp_in_arrays\tpct_of_array_bp\n")
        for c, bp in cls_bp.most_common():
            out.write(f"{args.sample}\t{c}\t{bp}\t{pct(bp)}\n")
    with open(os.path.join(args.outdir, "by_class_family.tsv"), "w") as out:
        out.write("sample\tclass_family\tbp_in_arrays\tpct_of_array_bp\n")
        for c, bp in cf_bp.most_common():
            out.write(f"{args.sample}\t{c}\t{bp}\t{pct(bp)}\n")
    with open(os.path.join(args.outdir, "by_family.tsv"), "w") as out:
        out.write("sample\tfamily\tclass_family\tlabel_source\tbp_in_arrays\tpct_of_array_bp\t"
                  "family_owned_bp\tfamily_share_in_arrays\n")
        for fam, bp in fam_bp.most_common(args.top):
            cf, origin = fam_label[fam]
            tot = total_fam.get(fam, 0)
            out.write(f"{args.sample}\t{fam}\t{cf}\t{origin or 'library'}\t{bp}\t{pct(bp)}\t{tot}\t"
                      f"{bp / tot:.3f}\n" if tot else f"{args.sample}\t{fam}\t{cf}\t{origin or 'library'}\t"
                      f"{bp}\t{pct(bp)}\t0\tNA\n")
    has_names = any(n != "." for ns in names.values() for n in ns)
    if has_names:
        per_motif = defaultdict(Counter)
        motif_len = Counter()
        for c, ivs in by_contig.items():
            for s, e, n in ivs:
                motif_len[n] += e - s
        idx_name = {}
        # clip_indexed's region index is per contig: rebuild names via clipped pieces
        for c, b, e, sc, (pl, i) in clipped:
            idx_name[(c, i)] = by_contig[c][i][2]
        for seg_bp, ((cf, _o, _f), i_region) in owned_segments(
                [(c, b, e, sc, (pl, (c, i))) for c, b, e, sc, (pl, i) in clipped]):
            per_motif[idx_name[i_region]][cf] += seg_bp
        with open(os.path.join(args.outdir, "by_motif.tsv"), "w") as out:
            out.write("sample\tmotif\tarray_bp\tmasked_bp\ttop_label\ttop_label_bp\tsatellite_bp\tunknown_bp\n")
            for n, bp in motif_len.most_common():
                labs = per_motif.get(n, Counter())
                top = labs.most_common(1)[0] if labs else ("unmasked", 0)
                sat = sum(v for k, v in labs.items() if collapse_class(k) == "Satellite")
                unk = sum(v for k, v in labs.items() if collapse_class(k) in ("Unknown", "Unknown_tandem"))
                out.write(f"{args.sample}\t{n}\t{bp}\t{sum(labs.values())}\t{top[0]}\t{top[1]}\t{sat}\t{unk}\n")

    lines = [f"sample\t{args.sample}",
             f"array_bp\t{array_bp} ({array_bp / 1e6:.1f} Mb, overlaps merged)",
             f"masked_by_repeatmasker\t{masked} ({pct(masked)}%)",
             f"unmasked\t{unmasked} ({pct(unmasked)}%)"]
    for c, bp in cls_bp.most_common():
        lines.append(f"class\t{c}\t{bp / 1e6:.2f} Mb\t{pct(bp)}%")
    lines.append("label source (curated table / applied cross-check / library + screens):")
    for (c, o), bp in sorted(origin_bp.items(), key=lambda kv: -kv[1])[:12]:
        lines.append(f"  {c}\t{o}\t{bp / 1e6:.2f} Mb")
    with open(os.path.join(args.outdir, "summary.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
