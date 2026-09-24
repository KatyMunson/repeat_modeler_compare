#!/usr/bin/env python3
"""Split one species' satellite-screen bp by where each hit motif came
from, using stage 02b's harmonized_summary.tsv (from the in-results
snapshot): the species' own species-coded motifs, motifs promoted to a
shared lowest-common-ancestor code (reported per lca_rank/lca_taxon),
other species' own motifs, or names absent from the summary (e.g. a
per-species override library). Never infers meaning from a name prefix
(MYX vs MYX2 can be genus or family depending on collision order).
bp are merged per contig within each origin group. Stdlib only."""

import argparse

from summarize_rm import merge_intervals


def load_summary(path):
    info = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            info[f[idx["promoted_name"]]] = (
                int(f[idx["n_species"]]), f[idx["species_codes"]], f[idx["lca_rank"]], f[idx["lca_taxon"]])
    return info


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-file", required=True, help="satellite-screen RepeatMasker .out")
    ap.add_argument("--harmonized-summary", required=True)
    ap.add_argument("--assembly-stats", required=True)
    ap.add_argument("--species", required=True, help="resolved species code")
    ap.add_argument("--species-id", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    info = load_summary(args.harmonized_summary)
    with open(args.assembly_stats) as fh:
        header = fh.readline().strip().split("\t")
        stats = dict(zip(header, fh.readline().strip().split("\t")))
    total_bp, non_n_bp = int(stats["total_bp"]), int(stats["non_n_bp"])

    groups = {}
    motifs = {}
    with open(args.out_file) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 11 or not f[0].isdigit() or not f[10].startswith("Satellite"):
                continue
            name = f[9]
            if name in info:
                n_species, codes, rank, taxon = info[name]
                if n_species > 1:
                    key = ("shared_promoted", rank, taxon)
                elif codes == args.species:
                    key = ("own_species", "species", args.species)
                else:
                    key = ("other_species", "species", codes)
            else:
                key = ("not_in_harmonized_summary", "NA", "NA")
            begin, end = sorted((int(f[5]), int(f[6])))
            groups.setdefault(key, {}).setdefault(f[4], []).append((begin, end))
            motifs.setdefault(key, set()).add(name)

    with open(args.out, "w") as out:
        out.write("species\tspecies_id\torigin\tlca_rank\tlca_taxon\tn_motifs_hit\tbp\tpct_total\tpct_non_n\n")
        for key in sorted(groups):
            bp = sum(merge_intervals(iv) for iv in groups[key].values())
            out.write(f"{args.species}\t{args.species_id}\t{key[0]}\t{key[1]}\t{key[2]}\t{len(motifs[key])}\t{bp}\t"
                      f"{100.0 * bp / total_bp:.4f}\t{100.0 * bp / non_n_bp:.4f}\n")


if __name__ == "__main__":
    main()
