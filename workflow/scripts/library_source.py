#!/usr/bin/env python3
"""Where does each species' masked sequence come from? Splits masked bp by
the library source of the family that holds it, from the combined
family_tandem.tsv (owned_bp: every base counted once, for its single
highest-scoring hit -- the rule class_composition.tsv uses):

  own_denovo        families discovered in this species (<species_id><sep>...)
  denovo:<species>  families discovered in another species of the run
  rm_builtin        RepeatMasker's own simple-repeat / low-complexity screen
  dfam              everything else (the Dfam export added to the library)

One row per species and source plus a `total` row per species:
  n_families  families holding any bp
  n_hits      .out hit lines (copies / fragments) of those families
  owned_bp, pct_of_masked, pct_non_n (bp / non-N assembly length from
  assembly_covariates.tsv, not a sum of per-family rounded percentages)
  bp_weighted_median_div  median of the families' median .out divergence,
                          weighted by owned bp (how old the matches are)

Counts follow the family that won each base at masking time. Dfam entries
are appended after the de novo families are clustered (never merged into
them), so where a de novo consensus fits better its bp count as de novo even
if Dfam holds the same element (dfam_overlap.tsv lists such families). A
family shared by both species is one cd-hit cluster named after its longest
member, so denovo:<other> includes shared families whose representative came
from the other species (library_membership.tsv).

With --dfam-matches (library/dfam_overlap/dfam_matches.tsv) each de novo
source is also split by the `dfam_match` column:
  all     every family of that source (the same numbers as without the flag)
  known   the de novo consensus matches a Dfam entry (cd-hit-est-2d,
          >= library.cdhit identity over coverage_short of the de novo family)
  novel   no such match
"known" is a lower bound on known material: diverged or partial matches
below the coverage cut count as novel, and a shared cluster counts as known
only if its representative (the family named in the library) matched.
Other rows carry `.` in dfam_match. Stdlib only."""

import argparse
import sys

BUILTIN_CLASSES = {"Simple_repeat", "Low_complexity"}


def source_of(family, cls, species, species_ids, sep):
    for sp in species_ids:
        if family.startswith(f"{sp}{sep}"):
            return "own_denovo" if sp == species else f"denovo:{sp}"
    if cls in BUILTIN_CLASSES:
        return "rm_builtin"
    return "dfam"


def weighted_median(pairs):
    pairs = sorted(pairs)
    total = sum(w for _v, w in pairs)
    acc = 0
    for v, w in pairs:
        acc += w
        if total and acc >= total / 2:
            return v
    return float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--family-tandem", required=True, help="combined family_tandem.tsv")
    ap.add_argument("--species-ids", nargs="+", required=True)
    ap.add_argument("--sep", default="_", help="library.species_prefix_sep")
    ap.add_argument("--assembly-covariates", required=True, help="assembly_covariates.tsv (species_id, non_n_bp)")
    ap.add_argument("--dfam-matches", help="dfam_matches.tsv from dfam_overlap (optional): split de novo "
                    "sources into known-in-Dfam and novel")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    non_n = {}
    with open(args.assembly_covariates) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            non_n[f[idx["species_id"]]] = int(float(f[idx["non_n_bp"]]))

    known = None
    if args.dfam_matches:
        known = set()
        with open(args.dfam_matches) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            i_fam = header.index("family")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) > i_fam:
                    known.add(f[i_fam].split("#", 1)[0])

    # (species, source, dfam_match) -> [n_families, n_hits, owned_bp, [(div, bp)]]
    agg = {}
    with open(args.family_tandem) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            owned = int(f[idx["owned_bp"]])
            if owned == 0:
                continue
            sp = f[idx["species"]]
            family = f[idx["family"]]
            src = source_of(family, f[idx["class"]], sp, args.species_ids, args.sep)
            keys = [(sp, src, "." if known is None or not src.startswith(("own_denovo", "denovo:")) else "all"),
                    (sp, "total", ".")]
            if known is not None and keys[0][2] == "all":
                keys.append((sp, src, "known" if family.split("#", 1)[0] in known else "novel"))
            for key in keys:
                a = agg.setdefault(key, [0, 0, 0, []])
                a[0] += 1
                a[1] += int(f[idx["n_hits"]])
                a[2] += owned
                try:
                    a[3].append((float(f[idx["median_div"]]), owned))
                except ValueError:
                    pass

    order = ["own_denovo"] + [f"denovo:{s}" for s in args.species_ids] + ["dfam", "rm_builtin", "total"]
    species = [s for s in args.species_ids if (s, "total", ".") in agg] + \
        sorted({s for s, _src, _m in agg} - set(args.species_ids))
    with open(args.out, "w") as out:
        out.write("species\tsource\tdfam_match\tn_families\tn_hits\towned_bp\tpct_of_masked\tpct_non_n\t"
                  "bp_weighted_median_div\n")
        for sp in species:
            masked = agg[(sp, "total", ".")][2]
            for src in order:
                for match in (".", "all", "known", "novel"):
                    if (sp, src, match) not in agg:
                        if match in ("known", "novel") and (sp, src, "all") in agg:
                            agg[(sp, src, match)] = [0, 0, 0, []]
                        else:
                            continue
                    n_fam, n_hits, bp, divs = agg[(sp, src, match)]
                    pct_nn = 100.0 * bp / non_n[sp] if non_n.get(sp) else float("nan")
                    out.write(f"{sp}\t{src}\t{match}\t{n_fam}\t{n_hits}\t{bp}\t"
                              f"{100.0 * bp / masked if masked else 0.0:.4f}\t"
                              f"{pct_nn:.4f}\t{weighted_median(divs):.1f}\n")
            own = agg.get((sp, "own_denovo", "all"), agg.get((sp, "own_denovo", "."), [0, 0, 0]))[2]
            dfam = agg.get((sp, "dfam", "."), [0, 0, 0])[2]
            msg = (f"[library_source] {sp}: own de novo {100.0 * own / masked:.1f}%, "
                   f"Dfam {100.0 * dfam / masked:.1f}% of {masked:,} masked bp")
            if known is not None:
                dn_known = sum(a[2] for (s_, src, m), a in agg.items() if s_ == sp and m == "known")
                dn_novel = sum(a[2] for (s_, src, m), a in agg.items() if s_ == sp and m == "novel")
                msg += (f"; de novo matching Dfam {100.0 * dn_known / masked:.1f}%, "
                        f"novel {100.0 * dn_novel / masked:.1f}%")
            print(msg, file=sys.stderr)


if __name__ == "__main__":
    main()
