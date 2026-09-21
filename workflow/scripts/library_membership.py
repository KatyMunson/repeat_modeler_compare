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
Stdlib only.
"""

import argparse
import re
import sys

MEMBER_RE = re.compile(r"^\d+\s+\d+(?:nt|aa), >(?P<name>.+)\.\.\.\s+(?P<rest>.*)$")


def species_of(name, sep):
    return name.split(sep, 1)[0]


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
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    clusters = parse_clstr(args.clstr)

    n_conflicts = 0
    with open(args.out, "w") as fh:
        fh.write(
            "cluster_id\trepresentative\trepresentative_class_family\tcategory\t"
            "n_members\tspecies_counts\tlabel_agreement\tdistinct_labels\n"
        )
        for cluster in clusters:
            members = cluster["members"]
            species_counts = {}
            labels = set()
            representative = None
            representative_label = None
            for m in members:
                sp = species_of(m["name"], args.sep)
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

            fh.write(
                f"{cluster['id']}\t{representative}\t{representative_label}\t{category}\t"
                f"{len(members)}\t{species_counts_str}\t{label_agreement}\t{distinct_labels_str}\n"
            )

    print(
        f"[library_membership] {len(clusters)} clusters, "
        f"{n_conflicts} with disagreeing Class/Family labels across members",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
