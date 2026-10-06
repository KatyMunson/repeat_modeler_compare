#!/usr/bin/env python3
"""Family counts through the discovery and library steps, per sample plus
an ALL row: how many families RepeatModeler's rounds and the LTR side
pipeline found, how many survived the RepeatModeler-style merge, how
RepeatClassifier labelled them, and how many ended up in clusters shared
with another sample versus clusters of their own after the cross-sample
cd-hit-est clustering (cluster_library).

Writes two tables:
  --out           one row per sample + ALL (wide)
  --by-class-out  the same family/cluster split per Class (the part of the
                  '#Class/Family' label before '/'), per sample + ALL

With --taxon (sample=taxon pairs from the manifest), clusters are also
counted as cross_taxon (members from more than one taxon). Two individuals
of one species sharing a family is then "shared" but not "cross_taxon", so
same-species samples don't inflate cross-species sharing.

Every sample's family count is checked against its members in the .clstr;
a mismatch means the inputs come from different runs and is an error.
Stdlib only.
"""

import argparse
import sys

from fasta_utils import iter_fasta, seq_id
from library_membership import class_family_of, parse_clstr, sample_of


def fasta_names(path):
    return [seq_id(header) for header, _ in iter_fasta(path)]


def class_of(name):
    return class_family_of(name).split("/", 1)[0]


def pct(part, whole):
    return f"{100.0 * part / whole:.2f}" if whole else "NA"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", nargs="+", required=True)
    ap.add_argument("--rounds-fa", nargs="+", required=True, help="per sample, same order as --sample")
    ap.add_argument("--ltr-fa", nargs="+", required=True)
    ap.add_argument("--merged-fa", nargs="+", required=True)
    ap.add_argument("--classified-fa", nargs="+", required=True)
    ap.add_argument("--clstr", required=True, help="cluster_library's shared_denovo.nr.fa.clstr")
    ap.add_argument("--membership", required=True, help="library_membership.tsv (for the Dfam count)")
    ap.add_argument("--sep", default="_")
    ap.add_argument("--taxon", nargs="*", default=[], help="sample=taxon (manifest); default: each sample its own")
    ap.add_argument("--out", required=True)
    ap.add_argument("--by-class-out", required=True)
    args = ap.parse_args()

    sample = args.sample
    for opt in ("rounds_fa", "ltr_fa", "merged_fa", "classified_fa"):
        if len(getattr(args, opt)) != len(sample):
            raise ValueError(f"--{opt.replace('_', '-')} needs one file per sample ({len(sample)})")
    known = set(sample)
    taxon = {sp: sp for sp in sample}
    taxon.update(dict(x.split("=", 1) for x in args.taxon))

    rows = {}
    for i, sp in enumerate(sample):
        merged_headers = [h for h, _ in iter_fasta(args.merged_fa[i])]
        classified = fasta_names(args.classified_fa[i])
        n_rounds = len(fasta_names(args.rounds_fa[i]))
        n_ltr = len(fasta_names(args.ltr_fa[i]))
        rows[sp] = {
            "rounds_families": n_rounds,
            "ltr_families": n_ltr,
            "merge_removed": n_rounds + n_ltr - len(merged_headers),
            "merged_families": len(merged_headers),
            "merged_ltr_families": sum(1 for h in merged_headers if seq_id(h).startswith("ltr-")),
            "putative_subfamilies": sum(1 for h in merged_headers if "putative subfamily of" in h),
            "classified_families": len(classified),
            "unknown_families": sum(1 for n in classified if class_of(n) == "Unknown"),
            "clusters": 0,
            "sample_only_clusters": 0,
            "shared_clusters": 0,
            "families_in_sample_only_clusters": 0,
            "families_in_shared_clusters": 0,
            "shared_clusters_label_conflict": 0,
            "cross_taxon_clusters": 0,
            "families_in_cross_taxon_clusters": 0,
        }

    by_class = {}  # (sample, class) -> [families, in_only, in_shared]
    totals = {"clusters": 0, "sample_only_clusters": 0, "shared_clusters": 0,
              "shared_clusters_label_conflict": 0, "cross_taxon_clusters": 0}
    for cluster in parse_clstr(args.clstr):
        members = [(sample_of(m["name"], args.sep, known), m["name"]) for m in cluster["members"]]
        present = {sp for sp, _ in members}
        shared = len(present) > 1
        cross = len({taxon[sp] for sp in present}) > 1
        totals["cross_taxon_clusters"] += cross
        conflict = shared and len({class_family_of(n) for _, n in members}) > 1
        totals["clusters"] += 1
        totals["shared_clusters" if shared else "sample_only_clusters"] += 1
        totals["shared_clusters_label_conflict"] += conflict
        for sp in present:
            row = rows[sp]
            row["clusters"] += 1
            row["shared_clusters" if shared else "sample_only_clusters"] += 1
            row["shared_clusters_label_conflict"] += conflict
            row["cross_taxon_clusters"] += cross
        for sp, name in members:
            rows[sp]["families_in_shared_clusters" if shared else "families_in_sample_only_clusters"] += 1
            rows[sp]["families_in_cross_taxon_clusters"] += cross
            for key in ((sp, class_of(name)), ("ALL", class_of(name))):
                counts = by_class.setdefault(key, [0, 0, 0])
                counts[0] += 1
                counts[2 if shared else 1] += 1

    for sp, row in rows.items():
        in_clusters = row["families_in_shared_clusters"] + row["families_in_sample_only_clusters"]
        if in_clusters != row["classified_families"]:
            raise ValueError(
                f"{sp}: {row['classified_families']} classified families but {in_clusters} in {args.clstr} "
                "-- inputs from different runs?"
            )

    # Unique names, so an older library_membership.tsv built from a
    # --add-reverse-complement export (each family twice) still counts families.
    dfam_names = set()
    with open(args.membership) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cat, rep = header.index("category"), header.index("representative")
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if fields[cat] == "dfam":
                dfam_names.add(fields[rep])
    n_dfam = len(dfam_names)

    all_row = {k: sum(r[k] for r in rows.values()) for k in next(iter(rows.values()))}
    all_row.update(totals)

    cols = ["rounds_families", "ltr_families", "merge_removed", "merged_families", "merged_ltr_families",
            "putative_subfamilies", "classified_families", "unknown_families", "clusters",
            "sample_only_clusters", "shared_clusters", "families_in_sample_only_clusters",
            "families_in_shared_clusters", "pct_families_in_shared_clusters",
            "shared_clusters_label_conflict", "cross_taxon_clusters", "families_in_cross_taxon_clusters",
            "pct_families_in_cross_taxon_clusters", "dfam_entries"]
    with open(args.out, "w") as fh:
        fh.write("sample\t" + "\t".join(cols) + "\n")
        for sp, row in list(rows.items()) + [("ALL", all_row)]:
            row = dict(row)
            row["pct_families_in_shared_clusters"] = pct(row["families_in_shared_clusters"],
                                                         row["classified_families"])
            row["pct_families_in_cross_taxon_clusters"] = pct(row["families_in_cross_taxon_clusters"],
                                                              row["classified_families"])
            row["dfam_entries"] = n_dfam if sp == "ALL" else "NA"
            fh.write(sp + "\t" + "\t".join(str(row[c]) for c in cols) + "\n")

    order = {sp: i for i, sp in enumerate(sample + ["ALL"])}
    with open(args.by_class_out, "w") as fh:
        fh.write("sample\tclass\tfamilies\tfamilies_in_sample_only_clusters\t"
                 "families_in_shared_clusters\tpct_in_shared_clusters\n")
        for (sp, cls), (n, n_only, n_shared) in sorted(by_class.items(), key=lambda kv: (order[kv[0][0]], kv[0][1])):
            fh.write(f"{sp}\t{cls}\t{n}\t{n_only}\t{n_shared}\t{pct(n_shared, n)}\n")

    a = all_row
    print(f"[discovery_summary] {a['classified_families']} families -> {a['clusters']} clusters "
          f"({a['shared_clusters']} shared, {a['sample_only_clusters']} sample-only); "
          f"{n_dfam} Dfam entries", file=sys.stderr)


if __name__ == "__main__":
    main()
