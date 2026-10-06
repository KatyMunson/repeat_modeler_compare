#!/usr/bin/env python3
"""What is in the sequence the LTR candidate tools skipped?

ltr_gather writes skipped_windows.tsv: the merged regions where an
LTRharvest or LTR_FINDER window (and its 50 kb salvage pieces) timed out.
Both tools slow down most on dense tandem arrays, so the working guess is
that skipped sequence is mostly satellite. This script checks that guess
against one RepeatMasker run (meant for the shared arm) and that run's
family_tandem.tsv.

Every base is given to its single highest-scoring .out hit, the rule
class_composition.tsv uses, so a base counts once. Rows per (sample, tool):

  level=total   name=skipped   bp = skipped bp, pct_of_skipped = 100
  level=total   name=masked    masked bp inside the skipped regions
  level=class   name=<class>   owned bp by class; Unknown families flagged
                               tandem_family are Unknown_tandem, as in
                               class_composition.tsv (curated and screen
                               labels are not applied here)
  level=class   name=unmasked  skipped bp no hit covers
  level=tandem  name=tandem_family / dispersed
                               owned bp of families family_tandem.py flags
                               tandem_family, whatever their class (an
                               LTR-labelled satellite counts as tandem)
  level=family  name=<family>  the --top-families families owning the most
                               skipped bp

pct_genome is the same quantity genome-wide (% of non-N bp), and
enrichment = pct_of_skipped / pct_genome. A skipped region that is mostly
tandem_family bp with high enrichment is an array the LTR tools could not
have found LTR candidates in anyway; a region that is mostly dispersed
LTR bp is discovery that was actually lost. Stdlib only."""

import argparse

from family_tandem import tandem_families
from intervals import clip
from summarize_rm import collapse_class, owned_segments, read_assembly_stats, relabel


def read_skipped(path):
    """{tool: {contig: sorted [(start, end)]}}, 0-based half-open."""
    skipped = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_tool, i_contig = header.index("tool"), header.index("contig")
        i_start, i_end = header.index("start"), header.index("end")
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            skipped.setdefault(f[i_tool], {}).setdefault(f[i_contig], []).append((int(f[i_start]), int(f[i_end])))
    for by_contig in skipped.values():
        for intervals in by_contig.values():
            intervals.sort()
    return skipped


def read_out_hits(path, tandem):
    """[(contig, begin, end, score, (class, family, is_tandem))], 1-based inclusive."""
    hits = []
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 11 or not f[0].isdigit():
                continue
            try:
                begin, end = sorted((int(f[5]), int(f[6])))
            except ValueError:
                continue
            family = f[9]
            cls = collapse_class(relabel(family, f[10], tandem)[0])
            hits.append((f[4], begin, end, int(f[0]), (cls, family, family in tandem)))
    return hits


def tally(segments):
    by_class, by_family, by_tandem = {}, {}, {}
    masked = 0
    for bp, (cls, family, is_tandem) in segments:
        masked += bp
        by_class[cls] = by_class.get(cls, 0) + bp
        by_family[family] = by_family.get(family, 0) + bp
        key = "tandem_family" if is_tandem else "dispersed"
        by_tandem[key] = by_tandem.get(key, 0) + bp
    return masked, by_class, by_family, by_tandem


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--skipped", required=True, help="ltr/skipped_windows.tsv (tool, contig, start, end)")
    ap.add_argument("--out-file", required=True, help="RepeatMasker .out (shared arm)")
    ap.add_argument("--tandem-table", required=True, help="family_tandem.tsv of the same run")
    ap.add_argument("--assembly-stats", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--top-families", type=int, default=10)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    _total_bp, non_n_bp = read_assembly_stats(args.assembly_stats)
    tandem = frozenset(tandem_families(args.tandem_table))
    hits = read_out_hits(args.out_file, tandem)
    skipped = read_skipped(args.skipped)

    g_masked, g_class, g_family, g_tandem = tally(owned_segments(hits))
    genome = {"total": {"masked": g_masked}, "class": g_class, "family": g_family, "tandem": g_tandem}
    genome["class"]["unmasked"] = non_n_bp - g_masked

    with open(args.out, "w") as out:
        out.write("sample\ttool\tlevel\tname\tbp\tpct_of_skipped\tpct_genome\tenrichment\n")

        def row(tool, level, name, bp, skipped_bp):
            pct = 100.0 * bp / skipped_bp if skipped_bp else 0.0
            g_bp = genome.get(level, {}).get(name)
            if level == "total" and name == "skipped":
                g_bp = None
            g_pct = 100.0 * g_bp / non_n_bp if g_bp is not None else None
            enr = pct / g_pct if g_pct else None
            out.write(f"{args.sample}\t{tool}\t{level}\t{name}\t{bp}\t{pct:.2f}\t"
                      f"{'NA' if g_pct is None else f'{g_pct:.4f}'}\t"
                      f"{'NA' if enr is None else f'{enr:.2f}'}\n")

        for tool in sorted(skipped):
            skipped_bp = sum(e - s for iv in skipped[tool].values() for s, e in iv)
            masked, by_class, by_family, by_tandem = tally(owned_segments(clip(hits, skipped[tool])))
            by_class["unmasked"] = skipped_bp - masked
            row(tool, "total", "skipped", skipped_bp, skipped_bp)
            row(tool, "total", "masked", masked, skipped_bp)
            for name in sorted(by_class, key=lambda k: -by_class[k]):
                row(tool, "class", name, by_class[name], skipped_bp)
            for name in ("tandem_family", "dispersed"):
                row(tool, "tandem", name, by_tandem.get(name, 0), skipped_bp)
            for name in sorted(by_family, key=lambda k: -by_family[k])[: args.top_families]:
                row(tool, "family", name, by_family[name], skipped_bp)


if __name__ == "__main__":
    main()
