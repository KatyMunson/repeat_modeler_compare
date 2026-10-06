#!/usr/bin/env python3
"""Per discovery round: how many of its families are NEW and how much of the
genome they own, as opposed to REFINEMENTS of families found earlier.

discovery_round_saturation.tsv counts the bp each round's families own. A
later round's sharper consensus of an earlier family takes over bases that
family already masked, so owned bp can stay high in a late round that adds
nothing. Here each round N is clustered with every earlier family (rounds
< N plus the LTR families) using RepeatModeler's own redundancy rule
(cd-hit-est -c 0.8 -aS 0.8 -g 1 -G 0, as in merge_families):
  refinement  a rnd-N family in a cluster with an earlier family
  new         every other rnd-N family (two rnd-N families clustering only
              with each other are both new)
bp come from the own-arm family_tandem.tsv (owned bp: every base given to
its best-scoring hit; a family missing from it, e.g. dropped as redundant by
merge_families or with no hits, owns 0).

Columns, per sample and round:
  n_families, n_new, n_refinement
  new_dispersed_pct_non_n   bp owned by new, non-tandem families: the number
                            that says whether another round still finds TE
                            families
  new_tandem_pct_non_n      bp owned by new tandem_family families
  refinement_pct_non_n      bp owned by refinements (not new coverage)
Stdlib only (cd-hit-est is external)."""

import argparse
import os
import re
import subprocess
import sys

from merge_families import parse_clusters

RND = re.compile(r"(?:^|_)(rnd-(\d+)_family-\d+)$")


def read_fasta(path):
    """[(id, header_line, seq)] in file order; id is the first header token."""
    recs, cur = [], None
    if not path or not os.path.exists(path):
        return recs
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                cur = [line[1:].split()[0], line, []]
                recs.append(cur)
            elif cur is not None:
                cur[2].append(line.strip())
    return [(i, h, "".join(s)) for i, h, s in recs]


def read_owned(path):
    """{rnd-N_family-M: (owned_bp, is_tandem)} from family_tandem.tsv."""
    out = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            m = RND.search(f[ix["family"]].split("#", 1)[0])
            if m:
                bp, tandem = out.get(m.group(1), (0, False))
                out[m.group(1)] = (bp + int(f[ix["owned_bp"]]), tandem or f[ix["tandem_family"]] == "True")
    return out


def refinements(cdhit, threads, workdir, earlier, current):
    """ids in `current` that cluster with any record of `earlier`."""
    if not earlier or not current:
        return set()
    fa = os.path.join(workdir, "combined.fa")
    with open(fa, "w") as out:
        for i, h, s in earlier + current:
            out.write(f">{i}\n{s}\n")
    cd_out = os.path.join(workdir, "cd-hit-out")
    cmd = [cdhit, "-aS", "0.8", "-c", "0.8", "-g", "1", "-G", "0", "-A", "80", "-M", "10000",
           "-d", "0", "-i", fa, "-o", cd_out, "-T", str(threads)]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    old = {i for i, _, _ in earlier}
    cur = {i for i, _, _ in current}
    refined = set()
    for members in parse_clusters(cd_out + ".clstr"):
        ids = {m[0] for m in members}
        if ids & old:
            refined |= ids & cur
    return refined


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rounds-fa", required=True, help="{sample}.rounds.consensi.fa")
    ap.add_argument("--ltr-fa", required=True, help="{sample}.ltrs.fa (may be empty)")
    ap.add_argument("--tandem-table", required=True, help="own-arm family_tandem.tsv")
    ap.add_argument("--assembly-stats", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--cdhit", default="cd-hit-est")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    with open(args.assembly_stats) as fh:
        header = fh.readline().strip().split("\t")
        non_n_bp = int(dict(zip(header, fh.readline().strip().split("\t")))["non_n_bp"])
    by_round = {}
    for rec in read_fasta(args.rounds_fa):
        m = RND.search(rec[0])
        if m:
            by_round.setdefault(int(m.group(2)), []).append(rec)
    ltr = [r for r in read_fasta(args.ltr_fa) if r[2]]
    owned = read_owned(args.tandem_table)
    os.makedirs(args.workdir, exist_ok=True)

    def pct(bp):
        return f"{100.0 * bp / non_n_bp:.4f}" if non_n_bp else "NA"

    with open(args.out, "w") as out:
        out.write("sample\tround\tn_families\tn_new\tn_refinement\tnew_dispersed_bp\tnew_tandem_bp\t"
                  "refinement_bp\tnew_dispersed_pct_non_n\tnew_tandem_pct_non_n\trefinement_pct_non_n\n")
        earlier = list(ltr)
        for n in sorted(by_round):
            current = by_round[n]
            refined = refinements(args.cdhit, args.threads, args.workdir, earlier, current)
            nd = nt = rb = 0
            for i, _, _ in current:
                bp, tandem = owned.get(i, (0, False))
                if i in refined:
                    rb += bp
                elif tandem:
                    nt += bp
                else:
                    nd += bp
            out.write(f"{args.sample}\trnd-{n}\t{len(current)}\t{len(current) - len(refined)}\t{len(refined)}\t"
                      f"{nd}\t{nt}\t{rb}\t{pct(nd)}\t{pct(nt)}\t{pct(rb)}\n")
            print(f"[round_novelty] {args.sample} rnd-{n}: {len(current)} families, {len(refined)} refinements "
                  f"of {len(earlier)} earlier; new dispersed {pct(nd)}% of non-N", file=sys.stderr)
            earlier += current


if __name__ == "__main__":
    main()
