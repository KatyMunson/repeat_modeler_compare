#!/usr/bin/env python3
"""Split a genome into N groups of WHOLE sequences, balanced by bp
(longest-first greedy: each sequence goes to the currently lightest
group). Sequences are never split across groups, so candidate
coordinates from different groups never need remapping and a group's
results are merged by plain concatenation. Writes one FASTA per group
(possibly empty when there are fewer sequences than groups) and a
manifest (group, contig, bp) so a slow or failing ltr_* group job can be
traced to its contigs. Stdlib only."""

import argparse
import heapq
import sys

from fasta_utils import iter_fasta, seq_id


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--outputs", nargs="+", required=True, help="one FASTA path per group")
    ap.add_argument("--manifest", required=True)
    args = ap.parse_args()

    lengths = []
    for header, seq in iter_fasta(args.fasta):
        lengths.append((len(seq), seq_id(header)))
    n = len(args.outputs)
    heap = [(0, g) for g in range(n)]
    assignment = {}
    for length, name in sorted(lengths, key=lambda t: (-t[0], t[1])):
        load, g = heapq.heappop(heap)
        assignment[name] = g
        heapq.heappush(heap, (load + length, g))

    outs = [open(p, "w") for p in args.outputs]
    for header, seq in iter_fasta(args.fasta):
        name = seq_id(header)
        out = outs[assignment[name]]
        out.write(f">{name}\n")
        for i in range(0, len(seq), 60):
            out.write(seq[i : i + 60] + "\n")
    for out in outs:
        out.close()

    with open(args.manifest, "w") as fh:
        fh.write("group\tcontig\tbp\n")
        for length, name in sorted(lengths, key=lambda t: (assignment[t[1]], -t[0])):
            fh.write(f"{assignment[name]}\t{name}\t{length}\n")
    loads = sorted(load for load, _ in heap)
    print(f"[group_genome] {len(lengths)} sequences into {n} groups; bp per group min={loads[0]} max={loads[-1]}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
