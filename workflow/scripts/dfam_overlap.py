#!/usr/bin/env python3
"""How many of each sample's de novo families are already-known Dfam
families: parse one cd-hit-est-2d .clstr per sample (db1 = the Dfam
export, db2 = that sample's prefixed families, same identity/coverage
thresholds as cluster_library) and report, per sample and Class, how
many families joined a Dfam cluster.

Writes:
  --out-summary  sample, class ("all" first), families, matching_dfam,
                 pct_matching_dfam, then matching_dfam split three ways:
                 class_agrees / class_differs (de novo family classified,
                 Dfam Class the same / different) and unknown_matched (de
                 novo family Unknown -- the match suggests a class);
                 per sample + ALL
  --out-matches  one row per matching de novo family: its label, the Dfam
                 family it clustered with, identity and strand

A match means the alignment covers >= --min-coverage of the DE NOVO
family: a fragment of a known element counts, a family that merely
contains a short Dfam entry (tRNA, MITE, solo LTR, satellite...) doesn't.
The rule's cd-hit-est-2d -s2 keeps those containment hits out of the
clustering itself; this check enforces the exact threshold and records the
coverage. Stdlib only.
"""

import argparse
import sys

from fasta_utils import iter_fasta, seq_id
from library_membership import class_family_of, parse_clstr


def class_of(name):
    return class_family_of(name).split("/", 1)[0]


def pct(part, whole):
    return f"{100.0 * part / whole:.2f}" if whole else "NA"


def parse_hit(rest, length):
    """'at 151:1650:1:1500/+/94.13%' (db2 start:end, db1 start:end) ->
    ('+', '94.13', fraction of the db2 sequence covered)."""
    fields = rest.split()[-1].split("/")
    if len(fields) != 3 or fields[0].count(":") != 3:
        raise ValueError(f"no alignment coordinates in .clstr hit {rest!r} (cd-hit-est-2d prints them by default)")
    a, b = (int(x) for x in fields[0].split(":")[:2])
    return fields[1], fields[2].rstrip("%"), (abs(b - a) + 1) / length


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", nargs="+", required=True)
    ap.add_argument("--families", nargs="+", required=True, help="prefixed FASTA per sample, same order")
    ap.add_argument("--clstr", nargs="+", required=True, help="cd-hit-est-2d .clstr per sample, same order")
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--out-matches", required=True)
    ap.add_argument("--min-coverage", type=float, default=0.8,
                    help="fraction of the de novo family the Dfam alignment must cover")
    args = ap.parse_args()
    if not len(args.sample) == len(args.families) == len(args.clstr):
        raise ValueError("--sample, --families and --clstr need one entry per sample")

    matches = []  # (sample, family, dfam_family, identity, strand, coverage)
    n_low_cov = 0
    counts = {}   # (sample, class) -> [families, matching, agrees, differs, unknown_matched]
    for sp, fa, clstr in zip(args.sample, args.families, args.clstr):
        names = [seq_id(h) for h, _ in iter_fasta(fa)]
        own = set(names)
        for name in names:
            for key in ((sp, "all"), (sp, class_of(name)), ("ALL", "all"), ("ALL", class_of(name))):
                counts.setdefault(key, [0, 0, 0, 0, 0])[0] += 1

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
                strand, identity, coverage = parse_hit(m["rest"], m["length"])
                if coverage < args.min_coverage:
                    n_low_cov += 1
                    continue
                matches.append((sp, m["name"], dfam[0]["name"], identity, strand, coverage))
                if class_of(m["name"]) == "Unknown":
                    col = 4
                else:
                    col = 2 if class_of(m["name"]) == class_of(dfam[0]["name"]) else 3
                for key in ((sp, "all"), (sp, class_of(m["name"])), ("ALL", "all"), ("ALL", class_of(m["name"]))):
                    counts[key][1] += 1
                    counts[key][col] += 1

    order = {sp: i for i, sp in enumerate(args.sample + ["ALL"])}
    with open(args.out_summary, "w") as fh:
        fh.write("sample\tclass\tfamilies\tmatching_dfam\tpct_matching_dfam\t"
                 "class_agrees\tclass_differs\tunknown_matched\n")
        for (sp, cls), (n, n_match, n_agree, n_differ, n_unknown) in sorted(
            counts.items(), key=lambda kv: (order[kv[0][0]], kv[0][1] != "all", kv[0][1])
        ):
            fh.write(f"{sp}\t{cls}\t{n}\t{n_match}\t{pct(n_match, n)}\t{n_agree}\t{n_differ}\t{n_unknown}\n")

    with open(args.out_matches, "w") as fh:
        fh.write("sample\tfamily\tclass_family\tdfam_family\tdfam_class_family\tidentity_pct\t"
                 "denovo_coverage_pct\tstrand\n")
        for sp, fam, dfam, identity, strand, coverage in matches:
            fh.write(f"{sp}\t{fam.split('#', 1)[0]}\t{class_family_of(fam)}\t"
                     f"{dfam.split('#', 1)[0]}\t{class_family_of(dfam)}\t{identity}\t"
                     f"{100 * coverage:.1f}\t{strand}\n")

    total = counts.get(("ALL", "all"), [0, 0, 0, 0, 0])
    print(f"[dfam_overlap] {total[1]} of {total[0]} de novo families match a Dfam family "
          f"({pct(total[1], total[0])}%); {n_low_cov} cd-hit hits dropped for covering "
          f"< {args.min_coverage:.0%} of the de novo family", file=sys.stderr)


if __name__ == "__main__":
    main()
