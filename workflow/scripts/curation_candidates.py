#!/usr/bin/env python3
"""Rank library families that deserve a manual recheck (profile, extend,
rebuild the consensus, group the pieces; README "Curating a family").

A family is a candidate when it holds >= --min-pct-masked of some species'
masked bp (owned_bp, combined family_tandem.tsv) AND at least one flag:

  unknown          still Unknown after the classify screens (unknown_tandem
                   when it is the tandem-carved Unknown)
  long_for_class   consensus longer than plausible for its label (LENGTH_MAX
                   below): usually a chimeric or composite consensus, or a
                   misclassification (a 4 kb "SINE", an 12 kb "Tc1")
  verify:<status>  class_disagreements.tsv status (disagree_order,
                   retroposon_vs_line, domain_conflict, rfam, ...). Maverick
                   domain_conflict is not flagged: Maverick integrases look
                   retroviral, so REXdb vs GyDB disagreeing is expected
  young            median .out divergence < --young-div: a recent burst,
                   biologically interesting and the easiest to curate
                   (near-identical copies)

Families already in an element group (--groups, summary.element_groups) are
listed with curated=<group> and ranked last, so finished work drops out.
Ranked by the largest pct_of_masked across species. Stdlib only."""

import argparse
import sys

# Longest-prefix match on class_family. Generous upper bounds for an
# autonomous element of each type; a consensus above them is suspect.
LENGTH_MAX = {
    "SINE": 1500,
    "Retroposon": 3000,
    "LINE": 9000,
    "LINE/Penelope": 6000,
    "PLE": 6000,
    "DNA": 12000,
    "DNA/TcMar": 5000,
    "DNA/hAT": 5000,
    "DNA/PiggyBac": 5000,
    "DNA/MULE": 6000,
    "DNA/PIF": 6000,
    "DNA/Maverick": 30000,
    "DNA/Polinton": 30000,
    "RC": 20000,
    "LTR": 15000,
    "LTR/DIRS": 8000,
    "LTR/Ngaro": 8000,
}
FLAG_STATUSES = {"disagree_order", "retroposon_vs_line", "domain_conflict", "rfam",
                 "rfam_other_rna", "domain_on_rna_label", "host_protein"}


def length_max(class_family):
    best = None
    for prefix, mx in LENGTH_MAX.items():
        # "DNA" covers "DNA/x"; "DNA/TcMar" covers "DNA/TcMar-Tc1"
        if class_family == prefix or class_family.startswith(prefix + "/") or \
                ("/" in prefix and class_family.startswith(prefix)):
            if best is None or len(prefix) > len(best[0]):
                best = (prefix, mx)
    return best


def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            yield dict(zip(header, f))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--family-tandem", required=True, help="combined family_tandem.tsv")
    ap.add_argument("--reclass", default="", help="unknown_reclassification.tsv")
    ap.add_argument("--disagreements", default="", help="class_disagreements.tsv")
    ap.add_argument("--groups", default="", help="element groups TSV (group, family)")
    ap.add_argument("--species-ids", nargs="*", default=[])
    ap.add_argument("--min-pct-masked", type=float, default=0.1)
    ap.add_argument("--young-div", type=float, default=3.0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    per = {}          # family -> species -> row
    masked = {}
    for r in read_tsv(args.family_tandem):
        sp = r["species"]
        masked[sp] = masked.get(sp, 0) + int(r["owned_bp"])
        per.setdefault(r["family"], {})[sp] = r
    species = [s for s in (args.species_ids or sorted(masked)) if s in masked]

    new_class, reclass_ev = {}, {}
    if args.reclass:
        for r in read_tsv(args.reclass):
            if r.get("new_class") and r["new_class"] != r.get("old_class"):
                new_class[r["family"]] = r["new_class"]
            reclass_ev[r["family"]] = r.get("evidence") or r.get("note") or ""
    status = {}
    if args.disagreements:
        for r in read_tsv(args.disagreements):
            status[r["family"]] = (r["status"], r.get("domain_call", ""))
    group = {}
    if args.groups:
        for r in read_tsv(args.groups):
            group[r["family"].split("#", 1)[0]] = r["group"]

    rows = []
    for fam, by_sp in per.items():
        any_row = next(iter(by_sp.values()))
        orig = any_row["class_family"]
        cls = new_class.get(fam, orig)
        pcts = {sp: 100.0 * int(by_sp[sp]["owned_bp"]) / masked[sp] if sp in by_sp else 0.0 for sp in species}
        top = max(pcts.values()) if pcts else 0.0
        if top < args.min_pct_masked:
            continue
        lead = max(by_sp, key=lambda sp: int(by_sp[sp]["owned_bp"]))
        cons_len = int(any_row["cons_len"]) if any_row["cons_len"].isdigit() else 0
        try:
            div = float(by_sp[lead]["median_div"])
        except ValueError:
            div = float("nan")
        tandem = any(r.get("tandem_family") == "True" for r in by_sp.values())

        flags = []
        if cls.split("/", 1)[0] == "Unknown":
            flags.append("unknown_tandem" if tandem else "unknown")
        lm = length_max(cls)
        if lm and cons_len > lm[1]:
            flags.append(f"long_for_class({cons_len}>{lm[1]})")
        st = status.get(fam)
        if st and st[0] in FLAG_STATUSES and not (st[0] == "domain_conflict" and "Maverick" in cls):
            flags.append(f"verify:{st[0]}" + (f"({st[1]})" if st[1] else ""))
        if div == div and div < args.young_div:
            flags.append(f"young({div:.1f}%)")
        if not flags:
            continue
        rows.append((fam in group, -top, fam, orig, cls, cons_len, flags, by_sp, pcts, div,
                     any_row.get("monomer_period", "NA"), group.get(fam, ""), reclass_ev.get(fam, "")))

    rows.sort()
    with open(args.out, "w") as out:
        out.write("rank\tfamily\trm_class_family\tclass_family\tcons_len\tflags\tcurated\t"
                  + "".join(f"owned_bp_{sp}\tpct_masked_{sp}\tn_hits_{sp}\tmedian_div_{sp}\t" for sp in species)
                  + "monomer_period\tclassify_evidence\n")
        for rank, (cur, _t, fam, orig, cls, cl, flags, by_sp, pcts, _d, mono, grp, ev) in enumerate(rows, 1):
            out.write(f"{rank}\t{fam}\t{orig}\t{cls}\t{cl}\t{','.join(flags)}\t{grp or '.'}\t")
            for sp in species:
                r = by_sp.get(sp)
                out.write(f"{r['owned_bp'] if r else 0}\t{pcts[sp]:.3f}\t{r['n_hits'] if r else 0}\t"
                          f"{r['median_div'] if r else 'NA'}\t")
            out.write(f"{mono}\t{ev or '.'}\n")
    n_open = sum(1 for r in rows if not r[0])
    print(f"[curation_candidates] {len(rows)} families flagged ({n_open} not yet curated) at "
          f">= {args.min_pct_masked}% of some species' masked bp", file=sys.stderr)


if __name__ == "__main__":
    main()
