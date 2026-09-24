#!/usr/bin/env python3
"""Gather per-group LTR candidate files into one whole-genome .scn for
LTRPipeline_from_scn, and report sequence skipped by window timeouts.

Candidates: LTRharvest-format lines from the patched LTR_HARVEST_parallel
and LTR_FINDER_parallel (-harvest_out). Both write 12 columns:
  s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) s(rLTR) e(rLTR) l(rLTR) sim seq-nr chr
where seq-nr is the index within THAT GROUP's FASTA, so indexes collide
across groups. Every line is rewritten to the 11-column format a
whole-genome `gt ltrharvest` (as RepeatModeler's LTRPipeline runs it)
produces, with seq-nr = 0-based index of `chr` in --genome, which is the
exact FASTA LTRPipeline_from_scn receives (it renames sequences seq1..seqN
in that order before LTR_retriever). The trailing name column is dropped
on purpose: LTR_retriever would otherwise try to match it against the
renamed seqN identifiers.

Validation (fails loudly): unknown sequence name; any coordinate outside
1..length; element start >= end; LTR boundaries outside the element.
Exact duplicate candidates (same seq, element and LTR coordinates) are
kept once.

Timeouts: --timeout-logs lines are `piece<TAB>bp<TAB>salvaged|skipped`.
Only `skipped` pieces lost sequence (salvaged ones were re-run as 50 kb
windows). Piece names are `<contig>_sub<N>` (window N of --window-size,
--overlap) or `<contig>_sub<N>_sub<M>` (salvage window M of 50 kb with
the salvage overlap). They are mapped back to genome intervals and
reported as a per-tool union (overlaps counted once).
Stdlib only."""

import argparse
import re
import sys

from fasta_utils import iter_fasta, merge_half_open, seq_id

SALVAGE_SIZE = 50000      # must match the vendored scripts' salvage call
SALVAGE_OVERLAP = 20000


def read_genome_index(path):
    order = {}
    lengths = {}
    for i, (header, seq) in enumerate(iter_fasta(path)):
        name = seq_id(header)
        order[name] = i
        lengths[name] = len(seq)
    return order, lengths


def parse_candidates(path, order, lengths, tool, seen, out_rows):
    n = 0
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split()
            if len(f) != 12:
                raise ValueError(f"{path}:{lineno}: expected 12 columns, got {len(f)}: {line.strip()!r}")
            name = f[11]
            if name not in order:
                raise ValueError(f"{path}:{lineno}: sequence {name!r} not in the genome FASTA")
            s, e, l, ls, le, ll, rs, re_, rl = (int(float(x)) for x in f[:9])
            length = lengths[name]
            for v in (s, e, ls, le, rs, re_):
                if not 1 <= v <= length:
                    raise ValueError(f"{path}:{lineno}: coordinate {v} outside 1..{length} on {name}")
            if not (s < e and s <= ls <= le <= rs <= re_ <= e):
                raise ValueError(f"{path}:{lineno}: inconsistent element/LTR coordinates: {line.strip()!r}")
            key = (name, s, e, ls, le, rs, re_)
            if key in seen:
                continue
            seen.add(key)
            out_rows.append((order[name], s, f"{s} {e} {l} {ls} {le} {ll} {rs} {re_} {rl} {f[9]} {order[name]}"))
            n += 1
    return n


def piece_interval(piece, lengths, window, overlap):
    """0-based half-open genome interval of a timed-out piece."""
    m = re.fullmatch(r"(.+)_sub(\d+)_sub(\d+)", piece)
    if m and m.group(1) in lengths:
        base, n1, n2 = m.group(1), int(m.group(2)), int(m.group(3))
        start = (n1 - 1) * (window - overlap) + (n2 - 1) * (SALVAGE_SIZE - SALVAGE_OVERLAP)
        return base, start, SALVAGE_SIZE
    m = re.fullmatch(r"(.+)_sub(\d+)", piece)
    if m and m.group(1) in lengths:
        return m.group(1), (int(m.group(2)) - 1) * (window - overlap), window
    raise ValueError(f"cannot map timed-out piece {piece!r} back to a genome sequence")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genome", required=True, help="the masked whole-genome FASTA given to LTRPipeline_from_scn")
    ap.add_argument("--harvest", nargs="*", default=[], help="LTR_HARVEST_parallel .scn per group")
    ap.add_argument("--finder", nargs="*", default=[], help="LTR_FINDER_parallel -harvest_out .scn per group")
    ap.add_argument("--timeout-logs", nargs="*", default=[], help="tool:path pairs, e.g. harvest:g0.timeouts.tsv")
    ap.add_argument("--window-size", type=int, required=True)
    ap.add_argument("--overlap", type=int, required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--out-scn", required=True)
    ap.add_argument("--out-skipped", required=True, help="tool, contig, start, end (0-based half-open, merged)")
    ap.add_argument("--out-summary", required=True, help="species, tool, n_candidates, skipped_bp")
    args = ap.parse_args()

    order, lengths = read_genome_index(args.genome)
    seen = set()
    rows = []
    counts = {"harvest": 0, "finder": 0}
    for path in args.harvest:
        counts["harvest"] += parse_candidates(path, order, lengths, "harvest", seen, rows)
    for path in args.finder:
        counts["finder"] += parse_candidates(path, order, lengths, "finder", seen, rows)
    rows.sort()

    with open(args.out_scn, "w") as fh:
        fh.write("# repeat_compare normalize_scn.py: LTR_HARVEST_parallel + LTR_FINDER_parallel candidates,\n")
        fh.write("# seq-nr = 0-based index in the masked genome FASTA\n")
        fh.write("# s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) s(rLTR) e(rLTR) l(rLTR) sim(LTRs) seq-nr\n")
        for _, _, line in rows:
            fh.write(line + "\n")

    skipped = {}
    for spec in args.timeout_logs:
        tool, path = spec.split(":", 1)
        with open(path) as fh:
            for line in fh:
                if not line.strip():
                    continue
                piece, _bp, action = line.rstrip("\n").split("\t")
                if action != "skipped":
                    continue
                contig, start, size = piece_interval(piece, lengths, args.window_size, args.overlap)
                end = min(start + size, lengths[contig])
                skipped.setdefault(tool, {}).setdefault(contig, []).append((start, end))

    skipped_bp = {}
    with open(args.out_skipped, "w") as fh:
        fh.write("tool\tcontig\tstart\tend\n")
        for tool in sorted(skipped):
            total = 0
            for contig in sorted(skipped[tool]):
                for s, e in merge_half_open(skipped[tool][contig]):
                    fh.write(f"{tool}\t{contig}\t{s}\t{e}\n")
                    total += e - s
            skipped_bp[tool] = total
    with open(args.out_summary, "w") as fh:
        fh.write("species\ttool\tn_candidates\tskipped_bp\n")
        for tool in ("harvest", "finder"):
            fh.write(f"{args.species}\t{tool}\t{counts[tool]}\t{skipped_bp.get(tool, 0)}\n")
    print(f"[normalize_scn] {len(rows)} unique candidates (harvest {counts['harvest']}, finder {counts['finder']}); "
          f"skipped bp {skipped_bp}", file=sys.stderr)


if __name__ == "__main__":
    main()
