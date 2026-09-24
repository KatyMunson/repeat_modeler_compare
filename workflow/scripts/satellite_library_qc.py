#!/usr/bin/env python3
"""Check every satellite-library motif against KM's operational
definition of a satellite, measured on this genome's satellite screen
(RepeatMasker -nolow -lib <satellite library>):

  length   monomer between --min-len and --max-len bp (75-2000: stage 02
           ran with min_period_length 75; 2000 is TRF's maximum period)
  copies   >= --min-copies genome-wide, estimated as merged hit bp /
           monomer length. (Stage 02's copy cutoffs are TRF per-ARRAY copy
           numbers -- copies inside one locus -- so they don't establish this.)
  tandem   >= --min-tandem-frac of the motif's bp in tandem arrays: hits of
           the same motif on the same contig and strand, chained when the
           gap is <= max(50 bp, 0.2 x monomer), whose span is >=
           --min-array-copies monomers. Mostly isolated hits = behaves like
           a dispersed (often TE-derived) repeat, whatever TRF saw.
  simple   <= --max-short-period-frac of the monomer covered by an exact
           period-1..10 self-repeat (a "satellite" that is really (CA)n
           counts as satellite here, because the screen runs with -nolow)

Report-only: nothing is filtered. Outputs one row per motif (--out) and a
one-row passing-only coverage summary (--passing-summary) so the headline
satellite % can be compared with the % from motifs that pass everything.
Only hits whose class starts with "Satellite" are used. Stdlib only."""

import argparse
import sys

from fasta_utils import iter_fasta
from summarize_rm import merge_intervals


def short_period_frac(seq, max_period=10, min_matches=8):
    """Fraction of seq covered by exact tandem self-repeats of period
    1..max_period: stretches where seq[i] == seq[i+p] for at least
    max(min_matches, 2p) consecutive positions (i.e. >= 3 exact copies of
    the unit; chance runs in random sequence are ~4^-8 or rarer)."""
    seq = seq.upper()
    n = len(seq)
    covered = bytearray(n)
    for p in range(1, max_period + 1):
        i = 0
        while i + p < n:
            if seq[i] != seq[i + p]:
                i += 1
                continue
            j = i
            while j + p < n and seq[j] == seq[j + p]:
                j += 1
            if (j - i) >= max(min_matches, 2 * p):  # positions i .. j+p-1 are periodic
                for k in range(i, j + p):
                    covered[k] = 1
            i = j + 1
    return sum(covered) / n if n else 0.0


def read_hits(path):
    """{motif: [(contig, strand, begin, end)]} for Satellite-class hits."""
    hits = {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 11 or not f[0].isdigit() or not f[10].startswith("Satellite"):
                continue
            begin, end = sorted((int(f[5]), int(f[6])))
            hits.setdefault(f[9], []).append((f[4], f[8], begin, end))
    return hits


def arrays_for(motif_hits, monomer_len):
    """Chain hits into arrays; returns [(contig, start, end, n_hits)]."""
    gap_max = max(50, int(0.2 * monomer_len))
    by_key = {}
    for contig, strand, b, e in motif_hits:
        by_key.setdefault((contig, strand), []).append((b, e))
    arrays = []
    for (contig, _strand), ivs in by_key.items():
        ivs.sort()
        cs, ce, k = ivs[0][0], ivs[0][1], 1
        for b, e in ivs[1:]:
            if b - ce - 1 <= gap_max:
                ce = max(ce, e)
                k += 1
            else:
                arrays.append((contig, cs, ce, k))
                cs, ce, k = b, e, 1
        arrays.append((contig, cs, ce, k))
    return arrays


def per_contig_bp(intervals):
    by_contig = {}
    for contig, b, e in intervals:
        by_contig.setdefault(contig, []).append((b, e))
    return sum(merge_intervals(iv) for iv in by_contig.values())


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--satellite-lib", required=True)
    ap.add_argument("--out-file", required=True, help="gathered satellite-screen RepeatMasker .out")
    ap.add_argument("--genomewide", required=True, help="satellite_genomewide.tsv (for the non-gap denominator)")
    ap.add_argument("--species", required=True)
    ap.add_argument("--species-id", required=True)
    ap.add_argument("--min-len", type=int, default=75)
    ap.add_argument("--max-len", type=int, default=2000)
    ap.add_argument("--min-copies", type=float, default=100)
    ap.add_argument("--min-array-copies", type=float, default=3)
    ap.add_argument("--min-tandem-frac", type=float, default=0.5)
    ap.add_argument("--max-short-period-frac", type=float, default=0.5)
    ap.add_argument("--out", required=True)
    ap.add_argument("--passing-summary", required=True)
    args = ap.parse_args()

    motifs = []
    for header, seq in iter_fasta(args.satellite_lib):
        name = header.split()[0].split("#", 1)[0]
        motifs.append((name, seq))
    hits = read_hits(args.out_file)
    unknown = sorted(set(hits) - {m for m, _ in motifs})
    if unknown:
        print(f"[satellite_library_qc] WARNING: {len(unknown)} hit names not in the library "
              f"(e.g. {unknown[:3]}) -- screen run with a different library?", file=sys.stderr)

    with open(args.genomewide) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        gw = dict(zip(header, fh.readline().rstrip("\n").split("\t")))
    nongap_bp = int(gw["nongap_bp"])

    passing_intervals = []
    n_pass = 0
    with open(args.out, "w") as out:
        out.write("species\tspecies_id\tmotif\tmonomer_len\tgenome_bp\test_copies\tn_hits\tn_arrays\t"
                  "n_isolated_hits\tlargest_array_copies\ttandem_bp\ttandem_frac\tshort_period_frac\t"
                  "pass_length\tpass_copies\tpass_tandem\tpass_not_simple\tpass_all\n")
        for name, seq in motifs:
            mlen = len(seq)
            mh = hits.get(name, [])
            genome_bp = per_contig_bp([(c, b, e) for c, _s, b, e in mh])
            arrays = arrays_for(mh, mlen) if mh else []
            tandem = [(c, b, e) for c, b, e, _k in arrays if (e - b + 1) / mlen >= args.min_array_copies]
            tandem_bp = min(per_contig_bp(tandem), genome_bp)
            est_copies = genome_bp / mlen if mlen else 0.0
            tandem_frac = tandem_bp / genome_bp if genome_bp else 0.0
            largest = max(((e - b + 1) / mlen for _c, b, e, _k in arrays), default=0.0)
            isolated = sum(1 for *_x, k in arrays if k == 1)
            spf = short_period_frac(seq)
            p_len = args.min_len <= mlen <= args.max_len
            p_copies = est_copies >= args.min_copies
            p_tandem = tandem_frac >= args.min_tandem_frac
            p_simple = spf <= args.max_short_period_frac
            p_all = p_len and p_copies and p_tandem and p_simple
            if p_all:
                n_pass += 1
                passing_intervals += [(c, b, e) for c, _s, b, e in mh]
            out.write(f"{args.species}\t{args.species_id}\t{name}\t{mlen}\t{genome_bp}\t{est_copies:.1f}\t"
                      f"{len(mh)}\t{len(arrays)}\t{isolated}\t{largest:.1f}\t{tandem_bp}\t{tandem_frac:.3f}\t"
                      f"{spf:.3f}\t{p_len}\t{p_copies}\t{p_tandem}\t{p_simple}\t{p_all}\n")

    passing_bp = per_contig_bp(passing_intervals)
    with open(args.passing_summary, "w") as out:
        out.write("species\tspecies_id\tn_motifs\tn_motifs_passing\tsatellite_bp_passing\t"
                  "satellite_pct_nongap_screen_passing\n")
        pct = 100.0 * passing_bp / nongap_bp if nongap_bp else 0.0
        out.write(f"{args.species}\t{args.species_id}\t{len(motifs)}\t{n_pass}\t{passing_bp}\t{pct:.4f}\n")
    print(f"[satellite_library_qc] {args.species}: {n_pass}/{len(motifs)} motifs pass all criteria; "
          f"passing motifs cover {pct:.2f}% of non-gap bp", file=sys.stderr)


if __name__ == "__main__":
    main()
