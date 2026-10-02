#!/usr/bin/env python3
"""Group library families that are pieces of one element (e.g. a
RepeatModeler internal-region family plus the LTR-pipeline families holding
its LTRs) and report the element's total masked footprint. Reporting only:
nothing is remasked, and each base still belongs to the one family that won
it (family_tandem.tsv owned_bp), so a group's bp is the sum of its members'.

Subcommands:
  members   from a blastn of a curated element consensus (query) against
            the library (subject), list library families that are pieces
            of it: >= --min-id identity over >= --min-cov of the library
            family's own length (merged HSPs). Writes a groups TSV:
              group  family  pident  family_cov  element_start  element_end
            blastn -query elem.fa -db library -outfmt \\
              "6 qseqid sseqid pident length qstart qend sstart send qlen slen evalue bitscore"
  report    sum owned_bp per species and group from the combined
            family_tandem.tsv. Writes one row per species x group plus one
            per member family (row_type = group / member).

Several groups can share one TSV (concatenate members outputs). A family
listed in two groups is counted in both and flagged on stderr.
Stdlib only."""

import argparse
import sys
from collections import defaultdict


def bare(name):
    return name.split("#", 1)[0]


def merged_len(ivs):
    total, cur_s, cur_e = 0, None, None
    for s, e in sorted(ivs):
        if cur_e is None or s > cur_e + 1:
            if cur_e is not None:
                total += cur_e - cur_s + 1
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        total += cur_e - cur_s + 1
    return total


def cmd_members(args):
    hsps = defaultdict(list)
    with open(args.blast) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12 or line.startswith("#"):
                continue
            pid, length = float(f[2]), int(f[3])
            qs, qe, ss, se, slen = int(f[4]), int(f[5]), int(f[6]), int(f[7]), int(f[9])
            hsps[bare(f[1])].append((pid, length, min(qs, qe), max(qs, qe), min(ss, se), max(ss, se), slen))
    rows = []
    for fam, hs in hsps.items():
        good = [h for h in hs if h[0] >= args.min_id]
        if not good:
            continue
        slen = good[0][6]
        cov = merged_len([(h[4], h[5]) for h in good]) / slen
        if cov < args.min_cov:
            continue
        pid = sum(h[0] * h[1] for h in good) / sum(h[1] for h in good)
        rows.append((fam, pid, cov, min(h[2] for h in good), max(h[3] for h in good)))
    rows.sort(key=lambda r: (r[3], -r[2]))
    with open(args.out, "w") as out:
        out.write("group\tfamily\tpident\tfamily_cov\telement_start\telement_end\n")
        for fam, pid, cov, a, b in rows:
            out.write(f"{args.group}\t{fam}\t{pid:.1f}\t{cov:.2f}\t{a}\t{b}\n")
    print(f"[family_groups] {args.group}: {len(rows)} library families >= {args.min_id}% identity over "
          f">= {args.min_cov:.0%} of their length", file=sys.stderr)


def cmd_report(args):
    members = defaultdict(list)
    seen = defaultdict(list)
    with open(args.groups) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ig, ifam = header.index("group"), header.index("family")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(ig, ifam) or not f[ifam]:
                continue
            members[f[ig]].append(bare(f[ifam]))
            seen[bare(f[ifam])].append(f[ig])
    for fam, gs in seen.items():
        if len(set(gs)) > 1:
            print(f"[family_groups] WARNING {fam} is in groups {sorted(set(gs))}; counted in each", file=sys.stderr)

    non_n = {}
    if args.assembly_covariates:
        with open(args.assembly_covariates) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            idx = {h: i for i, h in enumerate(header)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                non_n[f[idx["species_id"]]] = int(float(f[idx["non_n_bp"]]))

    fam_stats = defaultdict(dict)   # species -> family -> (owned_bp, n_hits, class_family)
    masked = defaultdict(int)
    with open(args.family_tandem) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            sp, owned = f[idx["species"]], int(f[idx["owned_bp"]])
            masked[sp] += owned
            fam_stats[sp][bare(f[idx["family"]])] = (owned, int(f[idx["n_hits"]]), f[idx["class_family"]])

    species = [s for s in (args.species_ids or sorted(fam_stats)) if s in fam_stats]
    with open(args.out, "w") as out:
        out.write("species\tgroup\trow_type\tfamily\tclass_family\tn_families\tn_hits\towned_bp\t"
                  "pct_of_masked\tpct_non_n\n")
        for sp in species:
            for group, fams in members.items():
                rows = [(fam,) + fam_stats[sp].get(fam, (0, 0, "")) for fam in dict.fromkeys(fams)]
                tot_bp = sum(r[1] for r in rows)
                tot_hits = sum(r[2] for r in rows)
                present = sum(1 for r in rows if r[1])

                def pct(bp, den):
                    return f"{100.0 * bp / den:.4f}" if den else "nan"
                out.write(f"{sp}\t{group}\tgroup\t.\t.\t{present}\t{tot_hits}\t{tot_bp}\t"
                          f"{pct(tot_bp, masked[sp])}\t{pct(tot_bp, non_n.get(sp))}\n")
                for fam, bp, nh, cls in sorted(rows, key=lambda r: -r[1]):
                    out.write(f"{sp}\t{group}\tmember\t{fam}\t{cls or '.'}\t{1 if bp else 0}\t{nh}\t{bp}\t"
                              f"{pct(bp, masked[sp])}\t{pct(bp, non_n.get(sp))}\n")
                print(f"[family_groups] {sp} {group}: {present}/{len(rows)} families, {tot_bp:,} bp "
                      f"({pct(tot_bp, masked[sp])}% of masked)", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("members")
    m.add_argument("--blast", required=True, help="blastn outfmt 6 (see above); query = element consensus")
    m.add_argument("--group", required=True, help="name for the element")
    m.add_argument("--min-id", type=float, default=80.0)
    m.add_argument("--min-cov", type=float, default=0.5, help="fraction of the library family's length")
    m.add_argument("--out", required=True)
    r = sub.add_parser("report")
    r.add_argument("--groups", required=True, help="members output (group, family columns)")
    r.add_argument("--family-tandem", required=True, help="combined family_tandem.tsv")
    r.add_argument("--assembly-covariates", default="")
    r.add_argument("--species-ids", nargs="*", default=[])
    r.add_argument("--out", required=True)
    args = ap.parse_args()
    {"members": cmd_members, "report": cmd_report}[args.cmd](args)


if __name__ == "__main__":
    main()
