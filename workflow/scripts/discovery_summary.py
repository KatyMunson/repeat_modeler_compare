#!/usr/bin/env python3
"""Family counts through the discovery and library steps, per species plus
an ALL row: how many families RepeatModeler's rounds and the LTR side
pipeline found, how many survived the RepeatModeler-style merge, how
RepeatClassifier labelled them, and how many ended up in clusters shared
with another species versus clusters of their own after the cross-species
cd-hit-est clustering (cluster_library).

Writes two tables:
  --out           one row per species + ALL (wide)
  --by-class-out  the same family/cluster split per Class (the part of the
                  '#Class/Family' label before '/'), per species + ALL

Every species' family count is checked against its members in the .clstr;
a mismatch means the inputs come from different runs and is an error.
Stdlib only.
"""

import argparse
import sys

from fasta_utils import iter_fasta, seq_id
from library_membership import class_family_of, parse_clstr, species_of


def fasta_names(path):
    return [seq_id(header) for header, _ in iter_fasta(path)]


def class_of(name):
    return class_family_of(name).split("/", 1)[0]


def pct(part, whole):
    return f"{100.0 * part / whole:.2f}" if whole else "NA"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--species", nargs="+", required=True)
    ap.add_argument("--rounds-fa", nargs="+", required=True, help="per species, same order as --species")
    ap.add_argument("--ltr-fa", nargs="+", required=True)
    ap.add_argument("--merged-fa", nargs="+", required=True)
    ap.add_argument("--classified-fa", nargs="+", required=True)
    ap.add_argument("--clstr", required=True, help="cluster_library's shared_denovo.nr.fa.clstr")
    ap.add_argument("--membership", required=True, help="library_membership.tsv (for the Dfam count)")
    ap.add_argument("--sep", default="_")
    ap.add_argument("--out", required=True)
    ap.add_argument("--by-class-out", required=True)
    args = ap.parse_args()

    species = args.species
    for opt in ("rounds_fa", "ltr_fa", "merged_fa", "classified_fa"):
        if len(getattr(args, opt)) != len(species):
            raise ValueError(f"--{opt.replace('_', '-')} needs one file per species ({len(species)})")
    known = set(species)

    rows = {}
    for i, sp in enumerate(species):
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
            "species_only_clusters": 0,
            "shared_clusters": 0,
            "families_in_species_only_clusters": 0,
            "families_in_shared_clusters": 0,
            "shared_clusters_label_conflict": 0,
        }

    by_class = {}  # (species, class) -> [families, in_only, in_shared]
    totals = {"clusters": 0, "species_only_clusters": 0, "shared_clusters": 0,
              "shared_clusters_label_conflict": 0}
    for cluster in parse_clstr(args.clstr):
        members = [(species_of(m["name"], args.sep, known), m["name"]) for m in cluster["members"]]
        present = {sp for sp, _ in members}
        shared = len(present) > 1
        conflict = shared and len({class_family_of(n) for _, n in members}) > 1
        totals["clusters"] += 1
        totals["shared_clusters" if shared else "species_only_clusters"] += 1
        totals["shared_clusters_label_conflict"] += conflict
        for sp in present:
            row = rows[sp]
            row["clusters"] += 1
            row["shared_clusters" if shared else "species_only_clusters"] += 1
            row["shared_clusters_label_conflict"] += conflict
        for sp, name in members:
            rows[sp]["families_in_shared_clusters" if shared else "families_in_species_only_clusters"] += 1
            for key in ((sp, class_of(name)), ("ALL", class_of(name))):
                counts = by_class.setdefault(key, [0, 0, 0])
                counts[0] += 1
                counts[2 if shared else 1] += 1

    for sp, row in rows.items():
        in_clusters = row["families_in_shared_clusters"] + row["families_in_species_only_clusters"]
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
            "species_only_clusters", "shared_clusters", "families_in_species_only_clusters",
            "families_in_shared_clusters", "pct_families_in_shared_clusters",
            "shared_clusters_label_conflict", "dfam_entries"]
    with open(args.out, "w") as fh:
        fh.write("species\t" + "\t".join(cols) + "\n")
        for sp, row in list(rows.items()) + [("ALL", all_row)]:
            row = dict(row)
            row["pct_families_in_shared_clusters"] = pct(row["families_in_shared_clusters"],
                                                         row["classified_families"])
            row["dfam_entries"] = n_dfam if sp == "ALL" else "NA"
            fh.write(sp + "\t" + "\t".join(str(row[c]) for c in cols) + "\n")

    order = {sp: i for i, sp in enumerate(species + ["ALL"])}
    with open(args.by_class_out, "w") as fh:
        fh.write("species\tclass\tfamilies\tfamilies_in_species_only_clusters\t"
                 "families_in_shared_clusters\tpct_in_shared_clusters\n")
        for (sp, cls), (n, n_only, n_shared) in sorted(by_class.items(), key=lambda kv: (order[kv[0][0]], kv[0][1])):
            fh.write(f"{sp}\t{cls}\t{n}\t{n_only}\t{n_shared}\t{pct(n_shared, n)}\n")

    a = all_row
    print(f"[discovery_summary] {a['classified_families']} families -> {a['clusters']} clusters "
          f"({a['shared_clusters']} shared, {a['species_only_clusters']} species-only); "
          f"{n_dfam} Dfam entries", file=sys.stderr)


if __name__ == "__main__":
    main()
