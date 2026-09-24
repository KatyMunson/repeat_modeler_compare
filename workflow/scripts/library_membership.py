#!/usr/bin/env python3
"""Parse a cd-hit-est .clstr file (from clustering all species' prefixed
RepeatModeler2 families together) into a per-cluster membership table:
representative, member count per species, a category (shared vs
{species}_only), and — since cd-hit-est keeps one arbitrary representative
per cluster and silently discards the rest — whether every member's
RepeatClassifier-assigned Class/Family label actually agrees. Two
independent per-species RepeatModeler runs can classify the same family
differently; a disagreement here means the final shared_library.fa entry's
class label is not necessarily consensus, which can bias
class_composition.tsv. Flagged, not resolved (consistent with the spec's
"warn, don't fail" treatment of caveats it doesn't want solved in code).

Each cluster row also carries `sources`: member counts by origin --
`rnd` (RepeatScout/RECON round family), `ltr` (LTR structural pipeline
family), `sat_relabeled` (a round family relabeled #Satellite by
satellite_relabel.py). Entries appended to the libraries after
clustering (Dfam export, satellite libraries) are listed as extra
one-member rows with category `dfam` / `satellite`, so the table
accounts for everything in shared_library.fa.
Stdlib only.
"""

import argparse
import re
import sys

MEMBER_RE = re.compile(r"^\d+\s+\d+(?:nt|aa), >(?P<name>.+)\.\.\.\s+(?P<rest>.*)$")


def species_of(name, sep, known_codes=None):
    code = name.split(sep, 1)[0]
    if known_codes and code not in known_codes:
        raise ValueError(f"family '{name}' has prefix '{code}', which is not one of the "
                         f"resolved species codes {sorted(known_codes)}")
    return code


def source_of(name):
    base = name.split("#", 1)[0]
    if "ltr-" in base and "_family-" in base:
        return "ltr"
    if class_family_of(name).startswith("Satellite"):
        return "sat_relabeled"
    return "rnd"


def appended_names(path):
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                yield line[1:].split()[0]


def class_family_of(name):
    if "#" in name:
        return name.split("#", 1)[1]
    return "Unknown"


def parse_clstr(path):
    clusters = []
    current = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">Cluster"):
                if current is not None:
                    clusters.append(current)
                current = {"id": line.split()[-1], "members": []}
                continue
            m = MEMBER_RE.match(line)
            if not m:
                raise ValueError(f"Unparseable .clstr line: {line!r}")
            name = m.group("name")
            is_rep = m.group("rest").strip() == "*"
            current["members"].append({"name": name, "is_rep": is_rep})
        if current is not None:
            clusters.append(current)
    return clusters


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--clstr", required=True)
    ap.add_argument("--sep", default="_", help="species_prefix_sep from config")
    ap.add_argument("--species-codes", nargs="*", default=[], help="resolved species codes (validates prefixes)")
    ap.add_argument("--dfam", help="Dfam export appended after clustering")
    ap.add_argument("--satellite-libs", nargs="*", default=[], help="satellite libraries appended after clustering")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    known = set(args.species_codes)

    clusters = parse_clstr(args.clstr)

    n_conflicts = 0
    with open(args.out, "w") as fh:
        fh.write(
            "cluster_id\trepresentative\trepresentative_class_family\tcategory\t"
            "n_members\tspecies_counts\tlabel_agreement\tdistinct_labels\tsources\n"
        )
        for cluster in clusters:
            members = cluster["members"]
            species_counts = {}
            source_counts = {}
            labels = set()
            representative = None
            representative_label = None
            for m in members:
                sp = species_of(m["name"], args.sep, known)
                src = source_of(m["name"])
                source_counts[src] = source_counts.get(src, 0) + 1
                species_counts[sp] = species_counts.get(sp, 0) + 1
                label = class_family_of(m["name"])
                labels.add(label)
                if m["is_rep"]:
                    representative = m["name"]
                    representative_label = label

            if representative is None:
                # cd-hit-est always marks exactly one member '*'; guard anyway.
                representative = members[0]["name"]
                representative_label = class_family_of(representative)

            category = "shared" if len(species_counts) > 1 else f"{next(iter(species_counts))}_only"
            label_agreement = len(labels) == 1
            if not label_agreement:
                n_conflicts += 1

            species_counts_str = ";".join(f"{sp}:{n}" for sp, n in sorted(species_counts.items()))
            distinct_labels_str = ";".join(sorted(labels))
            sources_str = ";".join(f"{k}:{v}" for k, v in sorted(source_counts.items()))

            fh.write(
                f"{cluster['id']}\t{representative}\t{representative_label}\t{category}\t"
                f"{len(members)}\t{species_counts_str}\t{label_agreement}\t{distinct_labels_str}\t{sources_str}\n"
            )

        appended = []
        if args.dfam:
            appended += [("dfam", n) for n in appended_names(args.dfam)]
        seen = set()
        for path in args.satellite_libs:
            for n in appended_names(path):
                if n not in seen:
                    seen.add(n)
                    appended.append(("satellite", n))
        for i, (category, name) in enumerate(appended, 1):
            label = class_family_of(name)
            fh.write(f"appended_{i}\t{name}\t{label}\t{category}\t1\tNA\tTrue\t{label}\t{category}:1\n")

    print(
        f"[library_membership] {len(clusters)} clusters, "
        f"{n_conflicts} with disagreeing Class/Family labels across members",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
