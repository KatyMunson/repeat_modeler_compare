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

Counts follow the family that won each base after the shared library was
clustered, so a Dfam entry merged into a de novo family counts as de novo
(dfam_overlap.tsv has the stricter view). Stdlib only."""

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
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    non_n = {}
    with open(args.assembly_covariates) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            non_n[f[idx["species_id"]]] = int(float(f[idx["non_n_bp"]]))

    # (species, source) -> [n_families, n_hits, owned_bp, [(div, bp)]]
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
            src = source_of(f[idx["family"]], f[idx["class"]], sp, args.species_ids, args.sep)
            for key in ((sp, src), (sp, "total")):
                a = agg.setdefault(key, [0, 0, 0, []])
                a[0] += 1
                a[1] += int(f[idx["n_hits"]])
                a[2] += owned
                try:
                    a[3].append((float(f[idx["median_div"]]), owned))
                except ValueError:
                    pass

    order = ["own_denovo"] + [f"denovo:{s}" for s in args.species_ids] + ["dfam", "rm_builtin", "total"]
    species = [s for s in args.species_ids if (s, "total") in agg] + \
        sorted({s for s, _ in agg} - set(args.species_ids))
    with open(args.out, "w") as out:
        out.write("species\tsource\tn_families\tn_hits\towned_bp\tpct_of_masked\tpct_non_n\t"
                  "bp_weighted_median_div\n")
        for sp in species:
            masked = agg[(sp, "total")][2]
            for src in order:
                if (sp, src) not in agg:
                    continue
                n_fam, n_hits, bp, divs = agg[(sp, src)]
                pct_nn = 100.0 * bp / non_n[sp] if non_n.get(sp) else float("nan")
                out.write(f"{sp}\t{src}\t{n_fam}\t{n_hits}\t{bp}\t{100.0 * bp / masked if masked else 0.0:.2f}\t"
                          f"{pct_nn:.4f}\t{weighted_median(divs):.1f}\n")
            own = agg.get((sp, "own_denovo"), [0, 0, 0])[2]
            dfam = agg.get((sp, "dfam"), [0, 0, 0])[2]
            print(f"[library_source] {sp}: own de novo {100.0 * own / masked:.1f}%, "
                  f"Dfam {100.0 * dfam / masked:.1f}% of {masked:,} masked bp", file=sys.stderr)


if __name__ == "__main__":
    main()
