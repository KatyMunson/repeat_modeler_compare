#!/usr/bin/env python3
"""Profile one repeat family: consensus structure plus what its genomic copies
look like. Diagnostic only, not run by the pipeline. Stdlib only.

  python3 workflow/scripts/family_profile.py --family Esto_rnd-1_family-332 \\
      --library results_v3/classify/library_consensi.fa \\
      --out-file results_v3/shared/Esto/repeatmasker/Esto.fa.out \\
      --genome results_v3/Esto/genome/Esto.fa --outdir retroposon_332

Writes to --outdir:
  summary.txt      consensus length/GC, 5'/3' tails (poly-A / poly-T), terminal
                   direct or inverted repeats (LTR / TIR), internal duplicated
                   segments (chimeric or tandem consensus), ORFs, and copy
                   statistics (insertions, full-length copies, divergence,
                   5'-truncation)
  orfs.faa         ORFs >= --min-orf-aa codons in all six frames (paste into
                   NCBI CD-search / blastp, or hmmscan against Pfam)
  coverage.tsv     per --bin bp of the consensus: insertions covering it and
                   their median divergence (a 5'-truncated retroelement is
                   flat-low then rises toward the 3' end)
  insertions.tsv   one row per insertion: .out fragments joined by RepeatMasker
                   ID, genomic span, consensus span, divergence
  contigs.tsv      per contig: insertions, bp, length, insertions per Mb
  full_length.fa   up to --max-copies near-full-length copies, oriented like
                   the consensus, with --flank bp each side (lowercase)
  tsd.tsv          per extracted copy: target-site duplication (longest exact
                   match, 4-30 bp, at the end of the left flank and the start of
                   the right flank, within --tsd-slop of the boundary) and the
                   3' poly-A run
"""

import argparse
import os
import statistics
import sys
from collections import defaultdict

from fasta_utils import iter_fasta, write_fasta

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
BASES = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: AA[16 * i + 4 * j + k]
         for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}


def revcomp(s):
    return s.translate(COMP)[::-1]


def translate(s):
    return "".join(CODON.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))


def orfs(seq, min_aa):
    """[(strand, frame, start_nt 1-based on +, end_nt, aa)] stop-to-stop ORFs."""
    out = []
    L = len(seq)
    for strand, s in (("+", seq), ("-", revcomp(seq))):
        for frame in range(3):
            prot = translate(s[frame:])
            pos = 0
            for piece in prot.split("*"):
                if len(piece) >= min_aa:
                    a = frame + 3 * pos
                    b = a + 3 * len(piece)
                    if strand == "+":
                        out.append((strand, frame, a + 1, b, piece))
                    else:
                        out.append((strand, frame, L - b + 1, L - a, piece))
                pos += len(piece) + 1
    return sorted(out, key=lambda o: -len(o[4]))


def tail_run(seq, base):
    n = 0
    for c in reversed(seq):
        if c == base:
            n += 1
        elif n and c == "N":
            continue
        else:
            break
    return n


def polya(core, right):
    """A run at the 3' end of a copy, continuing into the right flank (the
    .out boundary often stops inside the tail)."""
    n = tail_run(core, "A")
    for c in right:
        if c != "A":
            break
        n += 1
    return n


def a_rich_tail(seq, window=40):
    t = seq[-window:]
    return max((t.count(b) / len(t), b) for b in "ACGT") if t else (0.0, "-")


def longest_shared(a, b, k_min=8):
    """Longest exact substring of a also in b (simple k-mer seed + extend)."""
    best = (0, -1, -1)
    if min(len(a), len(b)) < k_min:
        return best
    index = defaultdict(list)
    for j in range(len(b) - k_min + 1):
        index[b[j:j + k_min]].append(j)
    for i in range(len(a) - k_min + 1):
        for j in index.get(a[i:i + k_min], ()):
            if i and j and a[i - 1] == b[j - 1]:
                continue  # not a maximal start
            n = k_min
            while i + n < len(a) and j + n < len(b) and a[i + n] == b[j + n]:
                n += 1
            if n > best[0]:
                best = (n, i, j)
    return best


def internal_duplications(seq, k=40):
    """Segments >= 60 bp present twice in the consensus (exact 40-mer seeds,
    either strand): a chimeric, tandemly duplicated or LTR-bearing consensus.
    Returns [(start, end, other_start, strand)] (1-based, first copy)."""
    L = len(seq)
    first = {}
    for i in range(L - k + 1):
        first.setdefault(seq[i:i + k], i)
    rc = revcomp(seq)
    runs = defaultdict(list)  # (strand, diagonal) -> positions of the first copy
    for i in range(L - k + 1):
        kmer = seq[i:i + k]
        j = first[kmer]
        if "N" not in kmer and j < i:
            runs[("+", i - j)].append(j)
    for i in range(L - k + 1):
        kmer = rc[i:i + k]
        j = first.get(kmer)
        pos = L - k - i  # where this kmer's reverse complement starts on +
        if j is not None and "N" not in kmer and j + k <= pos:
            runs[("-", j + pos)].append(j)
    out = []
    for (strand, diag), pos in runs.items():
        pos.sort()
        a = prev = pos[0]
        for p in pos[1:] + [None]:
            if p is not None and p == prev + 1:
                prev = p
                continue
            if prev + k - a >= 60:
                other = a + diag if strand == "+" else diag - prev
                out.append((a + 1, prev + k, other + 1, strand))
            if p is not None:
                a = prev = p
    return sorted(out)


def read_out(path, family):
    """Insertions of one family, .out fragments joined by the ID column.
    -> [dict(contig, start, end, strand, cons_start, cons_end, div, n_frag)]"""
    frags = defaultdict(list)
    with open(path) as fh:
        for n, line in enumerate(fh):
            f = line.split()
            if len(f) < 14 or not f[0].isdigit() or f[9] != family:
                continue
            b, e = sorted((int(f[5]), int(f[6])))
            if f[8] == "C":
                cs, ce = sorted((int(f[12]), int(f[13])))
            else:
                cs, ce = sorted((int(f[11]), int(f[12])))
            key = (f[4], f[14]) if len(f) > 14 else (f[4], f"line{n}")
            frags[key].append((b, e, f[8], cs, ce, float(f[1]), e - b + 1))
    ins = []
    for (contig, _id), fr in frags.items():
        bp = sum(x[6] for x in fr)
        ins.append(dict(
            contig=contig, start=min(x[0] for x in fr), end=max(x[1] for x in fr),
            strand="-" if fr[0][2] == "C" else "+",
            cons_start=min(x[3] for x in fr), cons_end=max(x[4] for x in fr),
            div=sum(x[5] * x[6] for x in fr) / bp, n_frag=len(fr)))
    return ins


def tsd(left, right, slop):
    """Longest exact match (4-30 bp) ending within `slop` of the end of `left`
    and starting within `slop` of the start of `right`."""
    best = ""
    for i in range(max(0, len(left) - slop - 30), len(left)):
        for j in range(0, slop + 1):
            n = 0
            while (i + n < len(left) and j + n < len(right) and n < 30
                   and left[i + n] == right[j + n] and left[i + n] != "N"):
                n += 1
            if n >= 4 and i + n >= len(left) - slop and n > len(best):
                best = left[i:i + n]
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--family", required=True)
    ap.add_argument("--library", required=True, help="FASTA holding the consensus (name or name#class)")
    ap.add_argument("--out-file", required=True, help="RepeatMasker .out")
    ap.add_argument("--genome", help="genome FASTA (for contig lengths, full-length copies, TSDs)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-orf-aa", type=int, default=100)
    ap.add_argument("--bin", type=int, default=100)
    ap.add_argument("--full-frac", type=float, default=0.9, help="consensus fraction a full-length copy spans")
    ap.add_argument("--max-copies", type=int, default=50)
    ap.add_argument("--flank", type=int, default=50)
    ap.add_argument("--tsd-slop", type=int, default=10)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    o = lambda name: os.path.join(args.outdir, name)

    cons = None
    for header, seq in iter_fasta(args.library):
        if header.split()[0].split("#", 1)[0] == args.family:
            cons, cons_header = seq.upper(), header
            break
    if cons is None:
        sys.exit(f"{args.family} not in {args.library}")
    L = len(cons)
    lines = [f"family\t{args.family}", f"header\t{cons_header}", f"consensus_len\t{L}",
             f"gc_pct\t{100.0 * (cons.count('G') + cons.count('C')) / max(1, L - cons.count('N')):.1f}"]

    # --- consensus structure
    frac3, b3 = a_rich_tail(cons)
    frac5, b5 = a_rich_tail(revcomp(cons))
    lines += [f"3prime_tail\tterminal runs A={tail_run(cons, 'A')} T={tail_run(cons, 'T')}; "
              f"last 40 bp {100 * frac3:.0f}% {b3}",
              f"5prime_tail\tfirst 40 bp (reverse complement) {100 * frac5:.0f}% {b5}"]
    end = min(500, L // 2)
    d_len, d_i, d_j = longest_shared(cons[:end], cons[-end:])
    i_len, i_i, i_j = longest_shared(cons[:end], revcomp(cons[-end:]))
    lines.append(f"terminal_direct_repeat\t{d_len} bp exact, shared by the first and last {end} bp "
                 "(LTR-like if >= ~100 and at the very ends)"
                 + (f" at 5' {d_i + 1} / 3' {L - end + d_j + 1}" if d_len else ""))
    lines.append(f"terminal_inverted_repeat\t{i_len} bp exact (TIR-like if >= ~10 at the very ends)"
                 + (f" at 5' {i_i + 1}" if i_len else ""))
    dups = internal_duplications(cons)
    lines.append(f"internal_duplications\t{len(dups)} segment(s) >= 60 bp present twice (40-mer exact)")
    for s, e, other, st in dups[:20]:
        lines.append(f"  dup\t{s}-{e} also at {other} ({st} strand)")
    found = orfs(cons, args.min_orf_aa)
    lines.append(f"orfs_ge_{args.min_orf_aa}aa\t{len(found)}")
    with open(o("orfs.faa"), "w") as fh:
        for k, (st, fr, a, b, aa) in enumerate(found, 1):
            lines.append(f"  orf{k}\t{st} strand, nt {a}-{b}, {len(aa)} aa")
            write_fasta(fh, f"{args.family}_orf{k} {st}:{a}-{b} {len(aa)}aa", aa)

    # --- genomic copies
    ins = read_out(args.out_file, args.family)
    if not ins:
        lines.append("insertions\t0 (family not in this .out)")
    else:
        spans = [x["cons_end"] - x["cons_start"] + 1 for x in ins]
        full = [x for x in ins if (x["cons_end"] - x["cons_start"] + 1) >= args.full_frac * L]
        divs = sorted(x["div"] for x in ins)
        reach5 = sum(1 for x in ins if x["cons_start"] <= 0.05 * L + 1)
        reach3 = sum(1 for x in ins if x["cons_end"] >= 0.95 * L)
        lines += [
            f"insertions\t{len(ins)} ({sum(x['n_frag'] for x in ins)} .out fragments)",
            f"genomic_bp\t{sum(x['end'] - x['start'] + 1 for x in ins)}",
            f"full_length\t{len(full)} spanning >= {args.full_frac:.0%} of the consensus",
            f"reach_5prime_end\t{reach5} ({100 * reach5 / len(ins):.1f}%) start within the first 5%",
            f"reach_3prime_end\t{reach3} ({100 * reach3 / len(ins):.1f}%) end within the last 5%",
            f"consensus_span_median\t{statistics.median(spans):.0f} bp",
            f"div_pct\tmedian {statistics.median(divs):.1f}; 10th pct {divs[len(divs) // 10]:.1f}; "
            f"90th pct {divs[9 * len(divs) // 10]:.1f}",
            f"div_full_length_median\t{statistics.median(x['div'] for x in full):.1f}" if full else
            "div_full_length_median\tNA",
        ]
        hist = defaultdict(int)
        for d in divs:
            hist[int(d)] += 1
        lines.append("div_histogram\t" + " ".join(f"{k}:{hist[k]}" for k in sorted(hist)))

        nb = (L + args.bin - 1) // args.bin
        cov = [[] for _ in range(nb)]
        for x in ins:
            for b in range((x["cons_start"] - 1) // args.bin, min(nb, (x["cons_end"] - 1) // args.bin + 1)):
                cov[b].append(x["div"])
        with open(o("coverage.tsv"), "w") as fh:
            fh.write("bin_start\tbin_end\tinsertions\tmedian_div\n")
            for b in range(nb):
                med = f"{statistics.median(cov[b]):.1f}" if cov[b] else "NA"
                fh.write(f"{b * args.bin + 1}\t{min(L, (b + 1) * args.bin)}\t{len(cov[b])}\t{med}\n")
        with open(o("insertions.tsv"), "w") as fh:
            fh.write("contig\tstart\tend\tstrand\tcons_start\tcons_end\tdiv\tn_frag\n")
            for x in sorted(ins, key=lambda x: (x["contig"], x["start"])):
                fh.write(f"{x['contig']}\t{x['start']}\t{x['end']}\t{x['strand']}\t{x['cons_start']}\t"
                         f"{x['cons_end']}\t{x['div']:.2f}\t{x['n_frag']}\n")

        per_contig = defaultdict(lambda: [0, 0])
        for x in ins:
            per_contig[x["contig"]][0] += 1
            per_contig[x["contig"]][1] += x["end"] - x["start"] + 1
        clen = {}
        picked = sorted(full, key=lambda x: x["div"])[:args.max_copies]
        want = defaultdict(list)
        for x in picked:
            want[x["contig"]].append(x)
        extracted = []
        if args.genome:
            for header, seq in iter_fasta(args.genome):
                name = header.split()[0]
                clen[name] = len(seq)
                for x in want.get(name, ()):
                    s0, e0 = x["start"] - 1, x["end"]
                    left = seq[max(0, s0 - args.flank):s0].upper()
                    core = seq[s0:e0].upper()
                    right = seq[e0:e0 + args.flank].upper()
                    if x["strand"] == "-":
                        left, core, right = revcomp(right), revcomp(core), revcomp(left)
                    extracted.append((x, left, core, right))
        with open(o("contigs.tsv"), "w") as fh:
            fh.write("contig\tlength\tinsertions\tbp\tinsertions_per_mb\n")
            for c, (n, bp) in sorted(per_contig.items(), key=lambda kv: -kv[1][0]):
                ln = clen.get(c)
                dens = f"{1e6 * n / ln:.2f}" if ln else "NA"
                fh.write(f"{c}\t{ln if ln else 'NA'}\t{n}\t{bp}\t{dens}\n")
        if clen:
            hit_contigs = len(per_contig)
            lines.append(f"contigs_with_insertions\t{hit_contigs} of {len(clen)}; "
                         f"genome-wide {1e6 * len(ins) / sum(clen.values()):.2f} insertions per Mb")
        if extracted:
            with open(o("full_length.fa"), "w") as fa, open(o("tsd.tsv"), "w") as ts:
                ts.write("copy\tcontig\tstart\tend\tstrand\tdiv\tcons_span\ttsd\ttsd_len\t3prime_polyA\n")
                n_tsd = 0
                for k, (x, left, core, right) in enumerate(extracted, 1):
                    name = f"{args.family}_copy{k} {x['contig']}:{x['start']}-{x['end']}({x['strand']}) div={x['div']:.1f}"
                    write_fasta(fa, name, left.lower() + core + right.lower())
                    t = tsd(left, right, args.tsd_slop)
                    n_tsd += len(t) >= 6
                    ts.write(f"copy{k}\t{x['contig']}\t{x['start']}\t{x['end']}\t{x['strand']}\t{x['div']:.2f}\t"
                             f"{x['cons_start']}-{x['cons_end']}\t{t}\t{len(t)}\t{polya(core, right)}\n")
                lines.append(f"tsd\t{n_tsd} of {len(extracted)} least-diverged full-length copies have a "
                             f">= 6 bp flanking direct repeat (tsd.tsv)")

    with open(o("summary.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
