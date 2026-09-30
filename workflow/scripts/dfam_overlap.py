#!/usr/bin/env python3
"""How many of each species' de novo families are already-known Dfam
families: parse one cd-hit-est-2d .clstr per species (db1 = the Dfam
export, db2 = that species' prefixed families, same identity/coverage
thresholds as cluster_library) and report, per species and Class, how
many families joined a Dfam cluster.

Writes:
  --out-summary  species, class ("all" first), families, matching_dfam,
                 pct_matching_dfam, class_agrees (matches whose de novo
                 Class equals the Dfam family's Class); per species + ALL
  --out-matches  one row per matching de novo family: its label, the Dfam
                 family it clustered with, identity and strand

A family carrying only a fragment of a Dfam family doesn't match (the
coverage threshold applies to the shorter sequence), so this counts
families that ARE a known family, not ones that contain a piece of one.
Stdlib only.
"""

import argparse
import sys

from fasta_utils import iter_fasta, seq_id
from library_membership import class_family_of, parse_clstr


def class_of(name):
    return class_family_of(name).split("/", 1)[0]


def pct(part, whole):
    return f"{100.0 * part / whole:.2f}" if whole else "NA"


def parse_hit(rest):
    """'at 301:3300:1:3000/+/90.27%' or 'at +/90.27%' -> ('+', '90.27')."""
    fields = rest.split()[-1].split("/")
    return fields[-2], fields[-1].rstrip("%")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--species", nargs="+", required=True)
    ap.add_argument("--families", nargs="+", required=True, help="prefixed FASTA per species, same order")
    ap.add_argument("--clstr", nargs="+", required=True, help="cd-hit-est-2d .clstr per species, same order")
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--out-matches", required=True)
    args = ap.parse_args()
    if not len(args.species) == len(args.families) == len(args.clstr):
        raise ValueError("--species, --families and --clstr need one entry per species")

    matches = []  # (species, family, dfam_family, identity, strand)
    counts = {}   # (species, class) -> [families, matching, class_agrees]
    for sp, fa, clstr in zip(args.species, args.families, args.clstr):
        names = [seq_id(h) for h, _ in iter_fasta(fa)]
        own = set(names)
        for name in names:
            for key in ((sp, "all"), (sp, class_of(name)), ("ALL", "all"), ("ALL", class_of(name))):
                counts.setdefault(key, [0, 0, 0])[0] += 1

        seen = set()
        for cluster in parse_clstr(clstr):
            dfam = [m for m in cluster["members"] if m["name"] not in own]
            hits = [m for m in cluster["members"] if m["name"] in own]
            if not hits:
                continue
            if len(dfam) != 1:
                raise ValueError(f"{clstr}: cluster {cluster['id']} has {len(dfam)} non-{sp} members, expected 1 Dfam family")
            for m in hits:
                if m["name"] in seen:
                    raise ValueError(f"{clstr}: {m['name']} is in more than one cluster")
                seen.add(m["name"])
                strand, identity = parse_hit(m["rest"])
                matches.append((sp, m["name"], dfam[0]["name"], identity, strand))
                agrees = class_of(m["name"]) == class_of(dfam[0]["name"])
                for key in ((sp, "all"), (sp, class_of(m["name"])), ("ALL", "all"), ("ALL", class_of(m["name"]))):
                    counts[key][1] += 1
                    counts[key][2] += agrees

    order = {sp: i for i, sp in enumerate(args.species + ["ALL"])}
    with open(args.out_summary, "w") as fh:
        fh.write("species\tclass\tfamilies\tmatching_dfam\tpct_matching_dfam\tclass_agrees\n")
        for (sp, cls), (n, n_match, n_agree) in sorted(
            counts.items(), key=lambda kv: (order[kv[0][0]], kv[0][1] != "all", kv[0][1])
        ):
            fh.write(f"{sp}\t{cls}\t{n}\t{n_match}\t{pct(n_match, n)}\t{n_agree}\n")

    with open(args.out_matches, "w") as fh:
        fh.write("species\tfamily\tclass_family\tdfam_family\tdfam_class_family\tidentity_pct\tstrand\n")
        for sp, fam, dfam, identity, strand in matches:
            fh.write(f"{sp}\t{fam.split('#', 1)[0]}\t{class_family_of(fam)}\t"
                     f"{dfam.split('#', 1)[0]}\t{class_family_of(dfam)}\t{identity}\t{strand}\n")

    total = counts.get(("ALL", "all"), [0, 0, 0])
    print(f"[dfam_overlap] {total[1]} of {total[0]} de novo families match a Dfam family "
          f"({pct(total[1], total[0])}%)", file=sys.stderr)


if __name__ == "__main__":
    main()
