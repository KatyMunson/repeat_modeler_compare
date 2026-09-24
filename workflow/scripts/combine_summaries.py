#!/usr/bin/env python3
"""Concatenate the per-(arm,species) chunks summarize_rm.py and
assembly_stats.py write into the final long-format comparison tables, and
compute arm_concordance.tsv (shared vs own pct + delta per species/class).
discovery_round_saturation.tsv / ltr_discovery.tsv: per-species chunks
concatenated. Stdlib only, no pandas.
"""

import argparse


def parse_manifest_tissue(manifest_path):
    tissue_by_species = {}
    with open(manifest_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                continue
            tissue_by_species[fields[0]] = fields[3]
    return tissue_by_species


def concat_chunks(chunk_paths, out_path):
    header = None
    with open(out_path, "w") as fh_out:
        for path in chunk_paths:
            with open(path) as fh_in:
                this_header = fh_in.readline()
                if header is None:
                    header = this_header
                    fh_out.write(header)
                elif this_header != header:
                    raise ValueError(f"Header mismatch in {path}: {this_header!r} != {header!r}")
                for line in fh_in:
                    fh_out.write(line)


def combine_assembly_covariates(chunk_paths, tissue_by_species, out_path):
    with open(out_path, "w") as fh_out:
        fh_out.write(
            "species_id\ttissue\tcontig_count\ttotal_bp\tn_bp\tnon_n_bp\tcontig_n50\tlargest_contig\n"
        )
        for path in chunk_paths:
            with open(path) as fh_in:
                fh_in.readline()  # header
                row = fh_in.readline().strip().split("\t")
            species_id = row[0]
            tissue = tissue_by_species.get(species_id, "unknown")
            fh_out.write("\t".join([species_id, tissue] + row[1:]) + "\n")


def build_arm_concordance(class_composition_path, out_path):
    # key: (species, class) -> {arm: (pct_total, pct_non_n)}
    by_key = {}
    with open(class_composition_path) as fh:
        header = fh.readline().strip().split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            arm = fields[idx["arm"]]
            species = fields[idx["species"]]
            cls = fields[idx["class"]]
            pct_total = float(fields[idx["pct_total"]])
            pct_non_n = float(fields[idx["pct_non_n"]])
            by_key.setdefault((species, cls), {})[arm] = (pct_total, pct_non_n)

    with open(out_path, "w") as fh_out:
        fh_out.write(
            "species\tclass\tpct_total_shared\tpct_total_own\tdelta_total\t"
            "pct_non_n_shared\tpct_non_n_own\tdelta_non_n\n"
        )
        for (species, cls), by_arm in sorted(by_key.items()):
            shared = by_arm.get("shared")
            own = by_arm.get("own")
            if shared is None or own is None:
                continue
            pct_total_shared, pct_non_n_shared = shared
            pct_total_own, pct_non_n_own = own
            fh_out.write(
                f"{species}\t{cls}\t{pct_total_shared:.4f}\t{pct_total_own:.4f}\t"
                f"{pct_total_shared - pct_total_own:.4f}\t{pct_non_n_shared:.4f}\t"
                f"{pct_non_n_own:.4f}\t{pct_non_n_shared - pct_non_n_own:.4f}\n"
            )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--class-chunks", nargs="+", required=True)
    ap.add_argument("--family-chunks", nargs="+", required=True)
    ap.add_argument("--divergence-chunks", nargs="+", required=True)
    ap.add_argument("--assembly-stats-chunks", nargs="+", required=True)
    ap.add_argument("--class-composition-out", required=True)
    ap.add_argument("--family-composition-out", required=True)
    ap.add_argument("--divergence-landscape-out", required=True)
    ap.add_argument("--assembly-covariates-out", required=True)
    ap.add_argument("--arm-concordance-out", required=True)
    ap.add_argument("--round-saturation-chunks", nargs="*", default=[])
    ap.add_argument("--round-saturation-out", required=True)
    ap.add_argument("--ltr-summary-chunks", nargs="*", default=[])
    ap.add_argument("--ltr-summary-out", required=True)
    args = ap.parse_args()

    tissue_by_species = parse_manifest_tissue(args.manifest)

    concat_chunks(args.class_chunks, args.class_composition_out)
    concat_chunks(args.family_chunks, args.family_composition_out)
    concat_chunks(args.divergence_chunks, args.divergence_landscape_out)
    combine_assembly_covariates(args.assembly_stats_chunks, tissue_by_species, args.assembly_covariates_out)
    build_arm_concordance(args.class_composition_out, args.arm_concordance_out)
    concat_chunks(args.round_saturation_chunks, args.round_saturation_out)
    concat_chunks(args.ltr_summary_chunks, args.ltr_summary_out)

    tissues = set(tissue_by_species.values())
    if "unknown" in tissues or len(tissues) > 1:
        print(
            f"[combine_summaries] WARNING: tissue values across species are "
            f"{sorted(tissues)} — see assembly_covariates.tsv / report before "
            f"comparing across species (germline vs soma genomes are different genomes).",
        )


if __name__ == "__main__":
    main()
