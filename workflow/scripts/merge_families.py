#!/usr/bin/env python3
"""Merge the RECON/RepeatScout round families with the LTR structural
families the way RepeatModeler 2.0.9 itself does after -LTRStruct
(port of the "Combine results between both pipelines" block in the
RepeatModeler script):

  1. combined.fa = ltr families + round consensi (same for .stk)
  2. cd-hit-est -aS 0.8 -c 0.8 -g 1 -G 0 -A 80 -M 10000 -T <threads>
  3. per cluster with >1 member:
       - any ltr-* member   -> every non-ltr member is redundant (dropped);
                               all ltr-* members are kept
       - rnd-* members only -> members other than the longest are kept
                               but tagged "[ putative subfamily of <longest> ]"
  4. write consensi.fa / families.stk without the redundant families,
     with the subfamily tags on the FASTA header and the Stockholm DE line.

Two deliberate differences from RepeatModeler 2.0.9, both edge cases:
  - RepeatModeler only evaluates a cluster when it reads the NEXT
    ">Cluster" line, so the last cluster in the .clstr is never evaluated.
    Here every cluster is.
  - RepeatModeler only rewrites consensi.fa when at least one family is
    redundant; if none is, its consensi.fa stays rounds-only and the LTR
    families are silently dropped. Here the combined set is always written.

The output is still unclassified; classify_families runs RepeatClassifier
on it, as RepeatModeler does next. Stdlib only (cd-hit-est is external).
"""

import argparse
import os
import re
import subprocess
import sys

CLSTR_MEMBER = re.compile(r"^\d+\s+(\d+)nt,\s*>((rnd|ltr)-\d+_family-\d+)")


def read_text(path):
    if path and os.path.exists(path):
        with open(path) as fh:
            return fh.read()
    return ""


def parse_clusters(clstr_path):
    clusters = []
    current = None
    with open(clstr_path) as fh:
        for line in fh:
            if line.startswith(">Cluster"):
                current = []
                clusters.append(current)
                continue
            m = CLSTR_MEMBER.match(line)
            if m:
                current.append((m.group(2), m.group(3), int(m.group(1))))
    return clusters


def resolve(clusters):
    redundant = set()
    subfamily_of = {}
    n_ltr = n_rnd = 0
    for members in clusters:
        n_ltr += sum(1 for _, t, _ in members if t == "ltr")
        n_rnd += sum(1 for _, t, _ in members if t == "rnd")
        if len(members) < 2:
            continue
        longest = {"ltr": (None, 0), "rnd": (None, 0)}
        for fam_id, fam_type, size in members:
            if size > longest[fam_type][1]:  # strict '>' as in RepeatModeler: first seen wins ties
                longest[fam_type] = (fam_id, size)
        if longest["ltr"][0]:
            redundant.update(fid for fid, t, _ in members if t != "ltr")
        elif longest["rnd"][0]:
            for fam_id, _, _ in members:
                if fam_id != longest["rnd"][0]:
                    subfamily_of[fam_id] = longest["rnd"][0]
    return redundant, subfamily_of, n_ltr, n_rnd


def write_fasta(combined_fa, redundant, subfamily_of, out_path):
    kept = 0
    with open(out_path, "w") as out:
        keep = False
        for line in combined_fa.splitlines(keepends=True):
            if line.startswith(">"):
                fam_id = line[1:].split()[0]
                keep = fam_id not in redundant
                if keep:
                    kept += 1
                    if fam_id in subfamily_of:
                        line = line.rstrip("\r\n") + f" [ putative subfamily of {subfamily_of[fam_id]} ]\n"
            if keep:
                out.write(line)
    return kept


def write_stk(combined_stk, redundant, subfamily_of, out_path):
    with open(out_path, "w") as out:
        fam_id = ""
        buf = []
        for line in combined_stk.splitlines(keepends=True):
            m = re.match(r"^#=GF\s+ID\s+(\S+)", line)
            if m:
                fam_id = m.group(1)
            if re.match(r"^#=GF\s+DE\s+\S", line) and fam_id in subfamily_of:
                line = line.rstrip("\r\n") + f" [ putative subfamily of {subfamily_of[fam_id]} ]\n"
            buf.append(line)
            if line.startswith("//"):
                if fam_id not in redundant:
                    out.write("".join(buf))
                buf = []
                fam_id = ""


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rounds-fa", required=True)
    ap.add_argument("--rounds-stk", required=True)
    ap.add_argument("--ltr-fa", required=True, help="may be empty (no LTR families found)")
    ap.add_argument("--ltr-stk", required=True)
    ap.add_argument("--cdhit", required=True, help="path to cd-hit-est")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out-fa", required=True)
    ap.add_argument("--out-stk", required=True)
    args = ap.parse_args()

    os.makedirs(args.workdir, exist_ok=True)
    combined_fa = read_text(args.ltr_fa) + read_text(args.rounds_fa)
    combined_stk = read_text(args.ltr_stk) + read_text(args.rounds_stk)
    combined_path = os.path.join(args.workdir, "combined.fa")
    with open(combined_path, "w") as fh:
        fh.write(combined_fa)

    if read_text(args.ltr_fa).strip():
        cd_out = os.path.join(args.workdir, "cd-hit-out")
        cmd = [args.cdhit, "-aS", "0.8", "-c", "0.8", "-g", "1", "-G", "0", "-A", "80", "-M", "10000",
               "-i", combined_path, "-o", cd_out, "-T", str(args.threads)]
        print("[merge_families] " + " ".join(cmd), file=sys.stderr)
        with open(os.path.join(args.workdir, "cd-hit-stdout"), "w") as log:
            subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)
        redundant, subfamily_of, n_ltr, n_rnd = resolve(parse_clusters(cd_out + ".clstr"))
    else:
        print("[merge_families] no LTR families -- rounds-only library", file=sys.stderr)
        redundant, subfamily_of, n_ltr, n_rnd = set(), {}, 0, combined_fa.count(">")

    kept = write_fasta(combined_fa, redundant, subfamily_of, args.out_fa)
    write_stk(combined_stk, redundant, subfamily_of, args.out_stk)
    print(f"[merge_families] {n_rnd} RepeatScout/RECON families, {n_ltr} LTRPipeline families, "
          f"removed {len(redundant)} redundant, {len(subfamily_of)} tagged putative subfamilies, "
          f"final family count = {kept}", file=sys.stderr)
    if kept == 0:
        raise ValueError("merge_families produced no families")


if __name__ == "__main__":
    main()
