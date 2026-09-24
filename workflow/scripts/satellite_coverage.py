#!/usr/bin/env python3
"""Satellite coverage from the satellite-only RepeatMasker screen (the
"upper bound": no competition from de novo/Dfam families).

Only hits whose Class starts with "Satellite" count (with -nolow the
screen emits no Simple_repeat/Low_complexity hits anyway; the filter
guards against a satellite_lib override carrying other classes). Hits
are merged per contig before summing, and N (gap) positions are
subtracted, so satellite bp is always non-gap bp.

Outputs:
  --bed        merged satellite intervals (0-based half-open), for ltr_mask
  --genomewide one row: species, species_id, tissue, total_bp, gap_bp,
               nongap_bp, satellite_bp, pct_total, pct_nongap
  --per-contig species, species_id, contig, contig_len, gap_bp,
               satellite_bp, pct_total, pct_nongap, genome_pct_nongap, flagged
               (flagged = pct_nongap >= --flag-pct; informational only)
Stdlib only."""

import argparse
import sys

from fasta_utils import intersect_bp, iter_fasta, merge_half_open, n_runs, seq_id
from summarize_rm import parse_out_file


def pct(num, den):
    return 100.0 * num / den if den else 0.0


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-file", required=True, help="gathered satellite-screen RepeatMasker .out")
    ap.add_argument("--fasta", required=True, help="the prepped genome the screen ran on")
    ap.add_argument("--species", required=True, help="resolved species code")
    ap.add_argument("--species-id", required=True)
    ap.add_argument("--tissue", required=True)
    ap.add_argument("--flag-pct", type=float, default=50.0)
    ap.add_argument("--bed", required=True)
    ap.add_argument("--genomewide", required=True)
    ap.add_argument("--per-contig", required=True)
    args = ap.parse_args()

    hits = {}
    n_hits = 0
    n_other = 0
    for query_seq, begin, end, class_family in parse_out_file(args.out_file):
        if not class_family.startswith("Satellite"):
            n_other += 1
            continue
        hits.setdefault(query_seq, []).append((begin - 1, end))  # .out is 1-based closed
        n_hits += 1
    if n_other:
        print(f"[satellite_coverage] ignored {n_other} non-Satellite hits", file=sys.stderr)

    rows = []
    seen = set()
    with open(args.bed, "w") as bed:
        for header, seq in iter_fasta(args.fasta):
            name = seq_id(header)
            seen.add(name)
            gaps = n_runs(seq)
            gap_bp = sum(e - s for s, e in gaps)
            sat = merge_half_open(hits.get(name, []))
            for s, e in sat:
                bed.write(f"{name}\t{s}\t{e}\n")
            sat_bp = sum(e - s for s, e in sat) - intersect_bp(sat, gaps)
            rows.append((name, len(seq), gap_bp, sat_bp))
    unknown = sorted(set(hits) - seen)
    if unknown:
        raise ValueError(
            f".out names not in the genome FASTA (screen ran on a different genome?): "
            f"{unknown[:5]} ({len(unknown)} total)"
        )

    total = sum(r[1] for r in rows)
    gap = sum(r[2] for r in rows)
    sat = sum(r[3] for r in rows)
    genome_pct_nongap = pct(sat, total - gap)
    with open(args.genomewide, "w") as fh:
        fh.write("species\tspecies_id\ttissue\ttotal_bp\tgap_bp\tnongap_bp\tsatellite_bp\tpct_total\tpct_nongap\n")
        fh.write(f"{args.species}\t{args.species_id}\t{args.tissue}\t{total}\t{gap}\t{total - gap}\t{sat}\t"
                 f"{pct(sat, total):.4f}\t{genome_pct_nongap:.4f}\n")
    n_flagged = 0
    with open(args.per_contig, "w") as fh:
        fh.write("species\tspecies_id\tcontig\tcontig_len\tgap_bp\tsatellite_bp\tpct_total\tpct_nongap\t"
                 "genome_pct_nongap\tflagged\n")
        for name, length, gap_bp, sat_bp in rows:
            p_nongap = pct(sat_bp, length - gap_bp)
            flagged = p_nongap >= args.flag_pct
            n_flagged += flagged
            fh.write(f"{args.species}\t{args.species_id}\t{name}\t{length}\t{gap_bp}\t{sat_bp}\t"
                     f"{pct(sat_bp, length):.4f}\t{p_nongap:.4f}\t{genome_pct_nongap:.4f}\t{flagged}\n")
    print(f"[satellite_coverage] {args.species}: {n_hits} hits, {sat} satellite bp "
          f"({genome_pct_nongap:.2f}% of non-gap), {n_flagged} contigs >= {args.flag_pct}%", file=sys.stderr)


if __name__ == "__main__":
    main()
