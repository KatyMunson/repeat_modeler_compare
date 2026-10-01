#!/usr/bin/env python3
"""Array-based tandem check for every library family in one RepeatMasker
run (meant for the shared arm): does each family sit in tandem arrays, like
a satellite, or is it dispersed, like a transposon?

Uses the same operational definition as the satellite QC that the removed
satellite arm ran (satellite_library_qc.py, tag satellite-arm-v1), so the
two can be compared motif-for-family:

  length   consensus between --min-len and --max-len bp (75-2000: TRF's
           min_period_length used upstream and TRF's maximum period)
  copies   >= --min-copies genome-wide, estimated as merged hit bp /
           consensus length
  tandem   >= --min-tandem-frac of the family's bp in tandem arrays: hits of
           the same family on the same contig and strand, chained when the
           gap is <= max(50 bp, 0.2 x consensus), whose span is >=
           --min-array-copies consensus lengths
  simple   <= --max-short-period-frac of the consensus covered by an exact
           period-1..10 self-repeat

`satellite_like` = all four pass. Fragments of one interrupted TE chain too,
but span about one consensus length, so they don't reach 3 copies. A
RepeatModeler consensus can be a multimer of the true monomer; chaining
works on hit positions, so that doesn't affect `tandem`, only `length`.

Per family, `owned_bp` is the bp the family holds once every base is given
to its single highest-scoring .out hit (the rule class_composition.tsv
uses), so owned_bp sums to the masked total. genome_bp, tandem_bp and the
array columns use all of the family's own hits, overlapped or not.
`median_div` is the bp-weighted median of .out's perc. div. column (raw
substitution %, not Kimura) and places a family on the landscape.

Writes one row per family (--out) and one row per class (--class-out):
class bp, and how much of it is held by satellite_like families.
Report-only: nothing is filtered or relabelled. Stdlib only."""

import argparse
import sys

from fasta_utils import iter_fasta
from summarize_rm import collapse_class, merge_intervals, owned_segments


def short_period_frac(seq, max_period=10, min_matches=8):
    """Fraction of seq covered by exact tandem self-repeats of period
    1..max_period (>= max(min_matches, 2p) consecutive periodic positions,
    i.e. >= 3 exact copies of the unit). Same as satellite_library_qc.py."""
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
            if (j - i) >= max(min_matches, 2 * p):
                for k in range(i, j + p):
                    covered[k] = 1
            i = j + 1
    return sum(covered) / n if n else 0.0


def _left(tok):
    return int(tok.strip("()"))


def read_out(path):
    """Returns (hits, class_of, cons_len_from_out).
    hits: list of (contig, begin, end, score, strand, family, pct_div).
    cons_len_from_out: per family, consensus length implied by the repeat
    begin/end/(left) columns (fallback when the library lacks the family)."""
    hits = []
    class_of = {}
    cons_len = {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 14 or not f[0].isdigit():
                continue
            try:
                begin, end = sorted((int(f[5]), int(f[6])))
                score = int(f[0])
                pct_div = float(f[1])
                if f[8] == "C":
                    implied = int(f[12]) + _left(f[11])
                else:
                    implied = int(f[12]) + _left(f[13])
            except ValueError:
                continue
            family = f[9]
            class_of.setdefault(family, f[10])
            if implied > cons_len.get(family, 0):
                cons_len[family] = implied
            hits.append((f[4], begin, end, score, f[8], family, pct_div))
    return hits, class_of, cons_len


def arrays_for(family_hits, cons_len):
    """Chain one family's hits into arrays: same contig and strand, gap <=
    max(50, 0.2 x consensus). Returns [(contig, start, end, n_hits)]."""
    gap_max = max(50, int(0.2 * cons_len))
    by_key = {}
    for contig, begin, end, strand in family_hits:
        by_key.setdefault((contig, strand), []).append((begin, end))
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


def weighted_median(pairs):
    """pairs: [(value, weight)]."""
    pairs = sorted(pairs)
    total = sum(w for _v, w in pairs)
    acc = 0
    for v, w in pairs:
        acc += w
        if acc >= total / 2:
            return v
    return float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-file", required=True, help="RepeatMasker .out")
    ap.add_argument("--library", required=True, help="library FASTA the .out was masked with")
    ap.add_argument("--assembly-stats", required=True)
    ap.add_argument("--arm", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--min-len", type=int, default=75)
    ap.add_argument("--max-len", type=int, default=2000)
    ap.add_argument("--min-copies", type=float, default=100)
    ap.add_argument("--min-array-copies", type=float, default=3)
    ap.add_argument("--min-tandem-frac", type=float, default=0.5)
    ap.add_argument("--max-short-period-frac", type=float, default=0.5)
    ap.add_argument("--major-min-copies", type=float, default=1000)
    ap.add_argument("--major-min-bp", type=int, default=100000)
    ap.add_argument("--out", required=True)
    ap.add_argument("--class-out", required=True)
    args = ap.parse_args()

    with open(args.assembly_stats) as fh:
        header = fh.readline().strip().split("\t")
        stats = dict(zip(header, fh.readline().strip().split("\t")))
    non_n_bp = int(stats["non_n_bp"])

    lib_seq = {}
    for header, seq in iter_fasta(args.library):
        lib_seq[header.split()[0].split("#", 1)[0]] = seq

    hits, class_of, out_cons_len = read_out(args.out_file)
    missing = sorted(set(class_of) - set(lib_seq))
    if missing:
        print(f"[family_tandem] WARNING: {len(missing)} .out families not in {args.library} "
              f"(e.g. {missing[:3]}); their consensus length comes from the .out columns "
              f"and short_period_frac is NA", file=sys.stderr)

    owned = {}
    for seg_bp, family in owned_segments((c, b, e, s, fam) for c, b, e, s, _st, fam, _d in hits):
        owned[family] = owned.get(family, 0) + seg_bp

    by_family = {}
    for contig, b, e, _s, strand, family, pct_div in hits:
        by_family.setdefault(family, []).append((contig, b, e, strand, pct_div))
    del hits

    class_rows = {}
    with open(args.out, "w") as out:
        out.write("arm\tspecies\tfamily\tclass_family\tclass\tcons_len\tn_hits\tgenome_bp\towned_bp\t"
                  "owned_pct_non_n\test_copies\tmedian_div\tn_arrays\tn_isolated_hits\tlargest_array_bp\t"
                  "largest_array_copies\ttandem_bp\ttandem_frac\tshort_period_frac\tpass_length\t"
                  "pass_copies\tpass_tandem\tpass_not_simple\tsatellite_like\tcopy_tier\n")
        for family in sorted(by_family, key=lambda f: -owned.get(f, 0)):
            fh_ = by_family[family]
            class_family = class_of[family]
            cls = collapse_class(class_family)
            seq = lib_seq.get(family)
            cons_len = len(seq) if seq else out_cons_len.get(family, 0)
            genome_bp = per_contig_bp([(c, b, e) for c, b, e, _st, _d in fh_])
            arrays = arrays_for([(c, b, e, st) for c, b, e, st, _d in fh_], cons_len) if cons_len else []
            tandem = [(c, b, e) for c, b, e, _k in arrays
                      if cons_len and (e - b + 1) / cons_len >= args.min_array_copies]
            tandem_bp = min(per_contig_bp(tandem), genome_bp)
            tandem_frac = tandem_bp / genome_bp if genome_bp else 0.0
            est_copies = genome_bp / cons_len if cons_len else 0.0
            largest_bp = max((e - b + 1 for _c, b, e, _k in arrays), default=0)
            largest_copies = largest_bp / cons_len if cons_len else 0.0
            n_isolated = sum(1 for *_x, k in arrays if k == 1)
            n_arrays = sum(1 for *_x, k in arrays if k > 1)
            med_div = weighted_median([(d, e - b + 1) for _c, b, e, _st, d in fh_])
            spf = short_period_frac(seq) if seq else None
            p_len = args.min_len <= cons_len <= args.max_len
            p_copies = est_copies >= args.min_copies
            p_tandem = tandem_frac >= args.min_tandem_frac
            p_simple = spf is None or spf <= args.max_short_period_frac
            sat_like = p_len and p_copies and p_tandem and p_simple
            if est_copies >= args.major_min_copies or genome_bp >= args.major_min_bp:
                tier = "major"
            elif p_copies:
                tier = "minor"
            else:
                tier = "below_floor"
            own = owned.get(family, 0)
            own_pct = 100.0 * own / non_n_bp if non_n_bp else 0.0
            spf_s = f"{spf:.3f}" if spf is not None else "NA"
            out.write(f"{args.arm}\t{args.species}\t{family}\t{class_family}\t{cls}\t{cons_len}\t{len(fh_)}\t"
                      f"{genome_bp}\t{own}\t{own_pct:.4f}\t{est_copies:.1f}\t{med_div:.1f}\t{n_arrays}\t"
                      f"{n_isolated}\t{largest_bp}\t{largest_copies:.1f}\t{tandem_bp}\t{tandem_frac:.3f}\t"
                      f"{spf_s}\t{p_len}\t{p_copies}\t{p_tandem}\t{p_simple}\t{sat_like}\t{tier}\n")

            row = class_rows.setdefault(cls, [0, 0, 0, 0, 0])
            row[0] += 1
            row[1] += own
            row[2] += tandem_bp
            if sat_like:
                row[3] += 1
                row[4] += own

    with open(args.class_out, "w") as out:
        out.write("arm\tspecies\tclass\tn_families\towned_bp\towned_pct_non_n\ttandem_bp_sum\t"
                  "n_satellite_like\tsatellite_like_owned_bp\tsatellite_like_pct_of_class\t"
                  "satellite_like_pct_non_n\n")
        for cls, (n_fam, own, tandem_sum, n_sat, sat_own) in sorted(class_rows.items()):
            out.write(f"{args.arm}\t{args.species}\t{cls}\t{n_fam}\t{own}\t"
                      f"{100.0 * own / non_n_bp if non_n_bp else 0.0:.4f}\t{tandem_sum}\t{n_sat}\t{sat_own}\t"
                      f"{100.0 * sat_own / own if own else 0.0:.2f}\t"
                      f"{100.0 * sat_own / non_n_bp if non_n_bp else 0.0:.4f}\n")

    n_sat = sum(r[3] for r in class_rows.values())
    sat_bp = sum(r[4] for r in class_rows.values())
    print(f"[family_tandem] {args.arm}/{args.species}: {n_sat} satellite-like families hold "
          f"{sat_bp:,} bp ({100.0 * sat_bp / non_n_bp if non_n_bp else 0.0:.2f}% of non-N)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
