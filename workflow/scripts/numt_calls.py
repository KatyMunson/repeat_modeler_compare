#!/usr/bin/env python3
"""Genome-wide NUMT calls from a blastn of the mitogenome against the
nuclear genome (README "Genome-wide NUMT calls"). Report-only.

Subcommands:

  double  write the mitogenome doubled (seq+seq, one record per input
          record, same name) so an alignment across the circular origin is
          one HSP; numt_blast then runs blastn -task dc-megablast with it as
          the query.
  call    blastn table -> NUMT calls for one sample.
  combine per-sample summary chunks -> summary/numt_summary.tsv.

`call`, in order:
  1. Fold query coordinates back modulo the mitogenome length L (qlen/2):
     mito_start/mito_end are 1-based on the mitogenome, read along the
     mitogenome's own strand; `wraps` = the hit crosses the origin.
  2. Mitochondrial contigs: hits at >= --mito-contig-id % covering
     >= --mito-contig-cov of a contig. Listed in mito_contigs.tsv, and no
     NUMT is called on them.
  3. Duplicates from the doubled query: the same genome locus is often
     reported once per query copy (whole, or as pieces at the copy edges).
     Hits are taken by decreasing bitscore, and one whose genome interval
     overlaps an already kept hit on the same contig and strand by
     >= --dup-overlap of its own length is dropped.
  4. Hit-level NUMTs: one per remaining hit (numt_hits.bed).
  5. Compound NUMTs (numts.bed): hits on the same contig and strand,
     neighbours in genome order, joined when the mitogenome coordinates
     continue (the next hit starts within --mito-merge-gap of where the
     previous one ended, gap or overlap, across the origin too) and the
     genome gap is <= --merge-gap. That is one insertion split by an indel
     or a later insertion. The largest gap's content is reported from the
     shared-arm .out (gap_fill: class:bp, best hit per base, the labels
     class_composition.tsv uses; `unmasked` = no hit).
     numt_gap_hist.tsv: the genome gap of every mito-colinear neighbour pair,
     whatever its size, binned, with what fills the gaps -- the evidence for
     choosing --merge-gap.
  6. Both levels keep calls with aligned bp >= --min-len (after merging,
     so a short piece can still join a compound call).

Coordinates are the prepped genome's (the names in the .out and every other
output); BED files are 0-based half-open with a '#' header line.

Sanity check (--satellite-calls, optional): families the satellite
cross-check proposes as Other/NUMT should sit inside called NUMTs. Their
.out bp in this sample and the part inside hit-level calls go to the log and
the summary. Stdlib only."""

import argparse
import bisect
import statistics
import sys

from fasta_utils import blast_subject_resolver, iter_fasta, merge_half_open, seq_id, write_fasta
from intervals import clip
from summarize_rm import CANONICAL_CLASSES, collapse_class, owned_segments, parse_out_file, read_assembly_stats

BLAST_FIELDS = "qseqid sseqid pident length qstart qend sstart send qlen slen evalue bitscore".split()
GAP_BINS = [0, 50, 100, 200, 500, 1000, 2000, 3000, 4000, 5000, 6000, 8000, 10000, 20000, 50000]
SUMMARY_CLASSES = sorted(CANONICAL_CLASSES | {"Other"}) + ["unmasked"]


# ---------------------------------------------------------------- double


def double(args):
    n = 0
    with open(args.out, "w") as out:
        for header, seq in iter_fasta(args.mito):
            write_fasta(out, seq_id(header), seq + seq)
            n += 1
            print(f"[numt] mitogenome {seq_id(header)}: {len(seq)} bp, doubled to {2 * len(seq)}")
    if not n:
        sys.exit(f"[numt] no sequence in {args.mito}")


# ---------------------------------------------------------------- hits


class Hit:
    __slots__ = ("mito", "contig", "start", "end", "strand", "pident", "aln_len", "mito_len",
                 "mito_start", "mito_end", "wraps", "bitscore")

    @property
    def bp(self):
        return self.end - self.start


def fold(qstart, qend, mito_len):
    """Doubled-query coordinates (1-based, qstart <= qend) -> (mito_start,
    mito_end, wraps) on the single mitogenome, 1-based. A hit longer than
    the mitogenome keeps wraps=True and its end folded too."""
    ms = (qstart - 1) % mito_len + 1
    me = (qend - 1) % mito_len + 1
    wraps = (qend - qstart + 1) > mito_len or me < ms
    return ms, me, wraps


def read_blast(path, resolve=lambda n: n):
    hits = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = dict(zip(BLAST_FIELDS, line.rstrip("\n").split("\t")))
            h = Hit()
            h.mito, h.contig = f["qseqid"], resolve(f["sseqid"])
            h.pident, h.aln_len, h.bitscore = float(f["pident"]), int(f["length"]), float(f["bitscore"])
            qlen = int(f["qlen"])
            if qlen % 2:
                sys.exit(f"[numt] {h.mito}: odd query length {qlen}; was the query doubled?")
            h.mito_len = qlen // 2
            qs, qe = int(f["qstart"]), int(f["qend"])
            h.mito_start, h.mito_end, h.wraps = fold(qs, qe, h.mito_len)
            ss, se = int(f["sstart"]), int(f["send"])
            h.strand = "+" if ss <= se else "-"
            h.start, h.end = min(ss, se) - 1, max(ss, se)  # 0-based half-open
            hits.append(h)
    return hits


def mito_contigs(hits, contig_len, min_cov, min_id):
    """{contig: (cov, identity)} for contigs covered >= min_cov by hits at
    >= min_id %. identity is weighted by alignment length."""
    by_contig = {}
    for h in hits:
        if h.pident >= min_id:
            by_contig.setdefault(h.contig, []).append(h)
    out = {}
    for contig, hs in by_contig.items():
        length = contig_len.get(contig)
        if not length:
            continue
        covered = sum(e - s for s, e in merge_half_open([(h.start, h.end) for h in hs]))
        cov = covered / length
        if cov >= min_cov:
            ident = sum(h.pident * h.aln_len for h in hs) / sum(h.aln_len for h in hs)
            out[contig] = (cov, ident)
    return out


def drop_duplicates(hits, min_overlap):
    """Best-bitscore first; drop a hit overlapping a kept one (same contig,
    strand) by >= min_overlap of its own genome length."""
    kept = {}
    out = []
    for h in sorted(hits, key=lambda x: (-x.bitscore, -x.bp, x.contig, x.start)):
        ivs = kept.setdefault((h.contig, h.strand), [])
        dup = False
        for s, e in ivs:
            if min(e, h.end) - max(s, h.start) >= min_overlap * h.bp:
                dup = True
                break
        if not dup:
            ivs.append((h.start, h.end))
            out.append(h)
    out.sort(key=lambda x: (x.contig, x.start, x.end, x.strand))
    return out


# ---------------------------------------------------------------- merging


def mito_step(prev, nxt):
    """Mito distance from where `prev` ends to where `nxt` starts, in the
    direction of the NUMT (genome order = prev, nxt): 0 = directly adjacent,
    > 0 a gap, < 0 an overlap. Taken modulo the mitogenome length into
    (-L/2, L/2], so a NUMT crossing the origin is continuous. On the minus
    strand, genome order runs backwards along the mitogenome."""
    L = prev.mito_len
    if prev.strand == "+":
        d = nxt.mito_start - prev.mito_end - 1
    else:
        d = prev.mito_start - nxt.mito_end - 1
    d %= L
    if d > L // 2:
        d -= L
    return d


def colinear(prev, nxt, mito_merge_gap):
    return (prev.mito == nxt.mito and prev.contig == nxt.contig and prev.strand == nxt.strand
            and abs(mito_step(prev, nxt)) <= mito_merge_gap)


def neighbour_pairs(hits, mito_merge_gap):
    """Mito-colinear neighbours: consecutive hits in genome order on the same
    contig and strand (opposite-strand hits in between don't break them)
    whose mito coordinates continue. Yields (prev, nxt, genome_gap)."""
    lanes = {}
    for h in hits:
        lanes.setdefault((h.contig, h.strand), []).append(h)
    for lane in lanes.values():
        lane.sort(key=lambda x: (x.start, x.end))
        for prev, nxt in zip(lane, lane[1:]):
            if colinear(prev, nxt, mito_merge_gap):
                yield prev, nxt, nxt.start - prev.end


def merge(hits, merge_gap, mito_merge_gap):
    """Chains of hits (lists, genome order). Each hit joins the chain whose
    last hit it continues (mito-colinear, genome gap <= merge_gap), choosing
    the closest such chain; otherwise it starts a new one."""
    lanes = {}
    for h in hits:
        lanes.setdefault((h.contig, h.strand), []).append(h)
    chains = []
    for lane in lanes.values():
        lane.sort(key=lambda x: (x.start, x.end))
        open_chains = []
        for h in lane:
            # lanes are sorted by start: a chain out of reach now stays so
            open_chains = [ch for ch in open_chains if h.start - ch[-1].end <= merge_gap]
            best = None
            for ch in open_chains:
                gap = h.start - ch[-1].end
                if colinear(ch[-1], h, mito_merge_gap) and (best is None or gap < best[0]):
                    best = (gap, ch)
            if best:
                best[1].append(h)
            else:
                chains.append([h])
                open_chains.append(chains[-1])
    return chains


class Compound:
    def __init__(self, chain):
        self.hits = chain
        first, last = chain[0], chain[-1]
        self.contig, self.strand, self.mito = first.contig, first.strand, first.mito
        self.start = min(h.start for h in chain)
        self.end = max(h.end for h in chain)
        self.aligned_bp = sum(h.bp for h in chain)
        self.identity = sum(h.pident * h.bp for h in chain) / self.aligned_bp
        # mito extent along the NUMT: genome order on +, reversed on -
        mfirst, mlast = (first, last) if self.strand == "+" else (last, first)
        self.mito_start, self.mito_end = mfirst.mito_start, mlast.mito_end
        self.wraps = any(h.wraps for h in chain) or self.mito_end < self.mito_start
        gaps = [(b.start - a.end, a.end, b.start) for a, b in zip(chain, chain[1:])]
        self.max_gap = max(gaps) if gaps else None


# ---------------------------------------------------------------- .out


class OutIndex:
    """Shared-arm .out hits with class labels, indexed per contig, for the
    class bp inside small regions (gaps) and the NUMTs themselves."""

    def __init__(self, path, tandem, reclass, curated):
        self.by_contig = {}
        for c, b, e, s, cf in parse_out_file(path, tandem, reclass, curated):
            self.by_contig.setdefault(c, []).append((b, e, s, collapse_class(cf[0])))
        self.starts, self.max_len = {}, {}
        for c, hs in self.by_contig.items():
            hs.sort()
            self.starts[c] = [h[0] for h in hs]
            self.max_len[c] = max(e - b + 1 for b, e, _s, _cls in hs)

    def _hits_in(self, contig, start, end):
        """.out hits overlapping [start, end) (0-based), as clip() input."""
        hs = self.by_contig.get(contig)
        if not hs:
            return
        starts = self.starts[contig]
        i = bisect.bisect_left(starts, start + 1 - self.max_len[contig])
        j = bisect.bisect_right(starts, end)
        for b, e, s, cls in hs[i:j]:
            if e > start:
                yield contig, b, e, s, cls

    def class_bp(self, regions):
        """regions: {contig: [(start, end)]} 0-based half-open, any order.
        Returns {class: bp} (best hit per base) + 'unmasked'."""
        out = {}
        total = masked = 0
        for contig, ivs in regions.items():
            for start, end in merge_half_open(ivs):
                total += end - start
                hits = list(self._hits_in(contig, start, end))
                for bp, cls in owned_segments(clip(hits, {contig: [(start, end)]})):
                    out[cls] = out.get(cls, 0) + bp
                    masked += bp
        out["unmasked"] = total - masked
        return out


def fill_string(class_bp):
    items = sorted(((c, bp) for c, bp in class_bp.items() if bp > 0), key=lambda x: -x[1])
    return ",".join(f"{c}:{bp}" for c, bp in items) or "."


def load_labels(args):
    tandem, reclass, curated = frozenset(), {}, {}
    if args.tandem_table:
        from family_tandem import tandem_families
        tandem = frozenset(tandem_families(args.tandem_table))
    if args.reclass_table:
        from reclassify_unknown import reclassified
        reclass = reclassified(args.reclass_table)
    if args.curated_table:
        from family_groups import curated_classes
        for path in args.curated_table:
            for fam, cf in curated_classes(path).items():
                curated.setdefault(fam, cf)
    return tandem, reclass, curated


def numt_families(path):
    fams = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_fam, i_prop = header.index("family"), header.index("proposed_class_family")
        i_conf = header.index("confidence")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > i_prop and f[i_prop] == "Other/NUMT":
                fams.append((f[i_fam], f[i_conf]))
    return fams


# ---------------------------------------------------------------- call


def read_fingerprint(path):
    lengths = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            k, v = line.rstrip("\n").split("\t")
            lengths[k] = int(v)
    return lengths


def quantiles(values):
    if not values:
        return ["NA"] * 5
    v = sorted(values)
    q = lambda p: v[min(len(v) - 1, int(round(p * (len(v) - 1))))]  # noqa: E731
    return [v[0], q(0.25), q(0.5), q(0.75), v[-1]]


def bin_label(gap):
    if gap <= 0:
        return 0
    i = bisect.bisect_left(GAP_BINS, gap)
    return i


def call(args):
    contig_len = read_fingerprint(args.fingerprint)
    _total_bp, non_n_bp = read_assembly_stats(args.assembly_stats)
    raw = read_blast(args.blast, blast_subject_resolver(contig_len))
    print(f"[numt] {args.sample}: {len(raw)} blastn HSPs")

    mito = mito_contigs(raw, contig_len, args.mito_contig_cov, args.mito_contig_id)
    with open(args.mito_contigs_out, "w") as out:
        out.write("sample\tcontig\tlength\tcov\tidentity\n")
        for c in sorted(mito):
            cov, ident = mito[c]
            out.write(f"{args.sample}\t{c}\t{contig_len[c]}\t{cov:.4f}\t{ident:.2f}\n")
            print(f"[numt] excluded as mitochondrial: {c} ({contig_len[c]} bp, cov {cov:.3f}, id {ident:.2f}%)")
    if not mito:
        print("[numt] no mitochondrial contig found (cov/identity thresholds)")

    nuclear = [h for h in raw if h.contig not in mito]
    hits = drop_duplicates(nuclear, args.dup_overlap)
    print(f"[numt] {len(nuclear)} HSPs on nuclear contigs, {len(hits)} after dropping doubled-query duplicates")

    tandem, reclass, curated = load_labels(args)
    out_index = OutIndex(args.out_file, tandem, reclass, curated)

    # numt_gap_hist.tsv: all mito-colinear neighbour pairs, any genome gap
    pairs = list(neighbour_pairs(hits, args.mito_merge_gap))
    hist = {}
    for prev, nxt, gap in pairs:
        b = bin_label(gap)
        row = hist.setdefault(b, {"n": 0, "single": 0, "fill": {}})
        row["n"] += 1
        if gap > 0:
            fill = out_index.class_bp({prev.contig: [(prev.end, nxt.start)]})
            top_cls, top_bp = max(fill.items(), key=lambda x: x[1])
            for c, bp in fill.items():
                row["fill"][c] = row["fill"].get(c, 0) + bp
            if top_cls != "unmasked" and top_bp >= 0.8 * gap:
                row["single"] += 1
    with open(args.gap_hist_out, "w") as out:
        out.write("sample\tgap_min\tgap_max\tn_pairs\tn_filled_one_class\tfill_bp\n")
        for b in range(len(GAP_BINS) + 1):
            lo = "<=0" if b == 0 else str(GAP_BINS[b - 1] + 1)
            hi = "0" if b == 0 else (str(GAP_BINS[b]) if b < len(GAP_BINS) else "inf")
            row = hist.get(b, {"n": 0, "single": 0, "fill": {}})
            out.write(f"{args.sample}\t{lo}\t{hi}\t{row['n']}\t{row['single']}\t{fill_string(row['fill'])}\n")
    print(f"[numt] {len(pairs)} mito-colinear neighbour pairs; "
          f"{sum(1 for p in pairs if p[2] <= args.merge_gap)} within merge_gap {args.merge_gap}")

    # hit level
    hit_calls = [h for h in hits if h.bp >= args.min_len]
    with open(args.hits_out, "w") as out:
        out.write("#contig\tstart\tend\tid\taligned_bp\tstrand\tidentity\tmito_start\tmito_end\tmito_wraps\n")
        for i, h in enumerate(hit_calls, 1):
            out.write(f"{h.contig}\t{h.start}\t{h.end}\t{args.sample}_numt_hit_{i}\t{h.bp}\t{h.strand}\t"
                      f"{h.pident:.2f}\t{h.mito_start}\t{h.mito_end}\t{'yes' if h.wraps else 'no'}\n")

    # compound
    comps = [Compound(ch) for ch in merge(hits, args.merge_gap, args.mito_merge_gap)]
    comps = sorted((c for c in comps if c.aligned_bp >= args.min_len), key=lambda c: (c.contig, c.start))
    with open(args.numts_out, "w") as out:
        out.write("#contig\tstart\tend\tid\taligned_bp\tstrand\tidentity\tmito_start\tmito_end\tmito_wraps\t"
                  "n_hits\tmax_gap\tgap_fill\n")
        for i, c in enumerate(comps, 1):
            if c.max_gap is None:
                gap, fill = ".", "."
            else:
                g, a_end, b_start = c.max_gap
                gap = str(g)
                fill = fill_string(out_index.class_bp({c.contig: [(a_end, b_start)]})) if g > 0 else "."
            out.write(f"{c.contig}\t{c.start}\t{c.end}\t{args.sample}_numt_{i}\t{c.aligned_bp}\t{c.strand}\t"
                      f"{c.identity:.2f}\t{c.mito_start}\t{c.mito_end}\t{'yes' if c.wraps else 'no'}\t"
                      f"{len(c.hits)}\t{gap}\t{fill}\n")
    print(f"[numt] {len(hit_calls)} hit-level NUMTs, {len(comps)} compound NUMTs "
          f"(merge_gap {args.merge_gap}, mito_merge_gap {args.mito_merge_gap}, min_len {args.min_len})")

    # summary chunk
    hit_regions = {}
    for h in hit_calls:
        hit_regions.setdefault(h.contig, []).append((h.start, h.end))
    numt_bp = sum(e - s for ivs in hit_regions.values() for s, e in merge_half_open(ivs))
    span_regions = {}
    for c in comps:
        span_regions.setdefault(c.contig, []).append((c.start, c.end))
    span_bp = sum(e - s for ivs in span_regions.values() for s, e in merge_half_open(ivs))
    overlap = out_index.class_bp(hit_regions)

    fam_cols = ["NA", "NA", "NA"]
    if args.satellite_calls:
        fams = numt_families(args.satellite_calls)
        fam_regions = {}
        for contig, b, e, fam in family_hits(args.out_file, {f for f, _ in fams}, mito):
            fam_regions.setdefault(fam, {}).setdefault(contig, []).append((b - 1, e))
        merged_hits = {c: merge_half_open(ivs) for c, ivs in hit_regions.items()}
        tot_bp = tot_in = 0
        print(f"[numt] sanity check: {len(fams)} families the cross-check calls Other/NUMT")
        for fam, conf in fams:
            regions = {c: merge_half_open(ivs) for c, ivs in fam_regions.get(fam, {}).items()}
            bp = sum(e - s for ivs in regions.values() for s, e in ivs)
            inside = sum(e - s for _c, s, e in intersect(regions, merged_hits))
            tot_bp += bp
            tot_in += inside
            frac = f"{inside / bp:.3f}" if bp else "NA"
            flag = "  <-- mostly outside called NUMTs" if bp and inside < 0.5 * bp else ""
            print(f"[numt]   {fam} ({conf}): {bp} bp in {args.sample} (nuclear contigs), "
                  f"{inside} inside hit-level NUMTs ({frac}){flag}")
        fam_cols = [str(len(fams)), str(tot_bp), str(tot_in)]
    else:
        print("[numt] sanity check skipped: no satellite cross-check calls (satellite_crosscheck not configured)")

    lens = [c.aligned_bp for c in comps]
    med_id = f"{statistics.median([c.identity for c in comps]):.2f}" if comps else "NA"
    with open(args.summary_out, "w") as out:
        out.write("\t".join(summary_header()) + "\n")
        row = [args.sample, "called", str(len(mito)), ",".join(sorted(mito)) or ".",
               str(len(hit_calls)), str(len(comps)),
               str(numt_bp), f"{100.0 * numt_bp / non_n_bp:.4f}" if non_n_bp else "NA",
               str(span_bp), med_id] + [str(x) for x in quantiles(lens)]
        row += [str(overlap.get(c, 0)) for c in SUMMARY_CLASSES]
        row += fam_cols
        out.write("\t".join(row) + "\n")


def family_hits(path, names, exclude_contigs):
    """(contig, begin, end, family) for .out hits of the named families
    (library name as in the .out), 1-based inclusive, off excluded contigs."""
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 11 or not f[0].isdigit() or f[9] not in names or f[4] in exclude_contigs:
                continue
            try:
                b, e = sorted((int(f[5]), int(f[6])))
            except ValueError:
                continue
            yield f[4], b, e, f[9]


def intersect(a, b):
    """Two {contig: sorted disjoint [(s, e)]} -> [(contig, s, e)] overlaps."""
    out = []
    for c, ivs in a.items():
        other = b.get(c, [])
        j = 0
        for s, e in ivs:
            while j < len(other) and other[j][1] <= s:
                j += 1
            k = j
            while k < len(other) and other[k][0] < e:
                lo, hi = max(s, other[k][0]), min(e, other[k][1])
                if lo < hi:
                    out.append((c, lo, hi))
                k += 1
    return out


def summary_header():
    return (["sample", "status", "n_mito_contigs", "mito_contigs", "n_numt_hits", "n_numts",
             "numt_bp", "numt_pct_non_n", "numt_span_bp", "median_identity",
             "len_min", "len_q25", "len_median", "len_q75", "len_max"]
            + [f"bp_overlap_{c}" for c in SUMMARY_CLASSES]
            + ["crosscheck_numt_families", "crosscheck_numt_family_bp", "crosscheck_numt_family_bp_in_numts"])


# ---------------------------------------------------------------- combine


def combine(args):
    header = summary_header()
    with open(args.out, "w") as out:
        out.write("\t".join(header) + "\n")
        for path in args.chunks:
            with open(path) as fh:
                fh.readline()
                for line in fh:
                    if line.strip():
                        out.write(line)
        for s in args.skipped:
            out.write("\t".join([s, "no_mitogenome"] + ["NA"] * (len(header) - 2)) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    d = sub.add_parser("double", help="write the mitogenome doubled (blastn query)")
    d.add_argument("--mito", required=True, help="MitoHiFi final_mitogenome.fasta")
    d.add_argument("--out", required=True)

    c = sub.add_parser("call", help="NUMT calls for one sample")
    c.add_argument("--blast", required=True, help="blastn -outfmt '6 " + " ".join(BLAST_FIELDS) + "'")
    c.add_argument("--sample", required=True)
    c.add_argument("--fingerprint", required=True, help="genome fingerprint.tsv (contig lengths)")
    c.add_argument("--assembly-stats", required=True)
    c.add_argument("--out-file", required=True, help="shared-arm RepeatMasker .out")
    c.add_argument("--tandem-table", help="family_tandem.tsv (Unknown_tandem, as summarize)")
    c.add_argument("--reclass-table", help="unknown_reclassification.tsv")
    c.add_argument("--curated-table", nargs="*", default=[], help="curated tables, first listing a family wins")
    c.add_argument("--satellite-calls", help="satellite_family_calls.tsv (sanity check; optional)")
    c.add_argument("--mito-contig-cov", type=float, default=0.8)
    c.add_argument("--mito-contig-id", type=float, default=98.0)
    c.add_argument("--dup-overlap", type=float, default=0.5)
    c.add_argument("--mito-merge-gap", type=int, default=100)
    c.add_argument("--merge-gap", type=int, default=500)
    c.add_argument("--min-len", type=int, default=50)
    c.add_argument("--hits-out", required=True)
    c.add_argument("--numts-out", required=True)
    c.add_argument("--gap-hist-out", required=True)
    c.add_argument("--mito-contigs-out", required=True)
    c.add_argument("--summary-out", required=True)

    m = sub.add_parser("combine", help="per-sample summaries -> numt_summary.tsv")
    m.add_argument("--chunks", nargs="*", default=[])
    m.add_argument("--skipped", nargs="*", default=[], help="samples without a mitogenome")
    m.add_argument("--out", required=True)

    args = ap.parse_args()
    {"double": double, "call": call, "combine": combine}[args.cmd](args)


if __name__ == "__main__":
    main()
