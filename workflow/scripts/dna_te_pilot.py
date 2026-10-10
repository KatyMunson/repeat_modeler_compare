#!/usr/bin/env python3
"""Structural TIR / Helitron pilot (README "Structural TIR / Helitron
pilot"). Report-only: the library is not changed.

  gather   per-group candidates from dna_te_windows.py -> one FASTA per type
           (candidates reported twice by overlapping windows dropped),
           candidates.tsv, and skipped_windows.tsv (timed-out or failed
           windows, merged per tool and reason, as ltr/skipped_windows.tsv).
  cluster  per sample and type: cd-hit-est the candidates (library.cdhit
           settings, both strands); cd-hit-est-2d the cluster
           representatives against the sample's own -families.fa and
           against the shared library; split the unmatched
           representatives into length-balanced chunks.
  blast_chunk  blastn (dc-megablast) one chunk against the genome; reduce
           the HSPs on disk to per-query copies and merged intervals.
  finalize per sample: chunk summaries -> one row per type (pilot.tsv) and
           the per-cluster table with genome copies.
  combine  per-sample rows -> summary/dna_te_pilot.tsv, with the decision
           rule on its first line.

The bp the unmatched clusters would mask is a blastn proxy: the union of
their dc-megablast HSPs on the genome, which is less sensitive than
RepeatMasker -s (old, diverged copies are missed), so it is a lower bound.
`new_bp` is the part of it the shared-arm .out does not mask already; the
decision uses new_bp. `*_repeated` counts only representatives with
>= --min-copies genomic loci besides their own (informational for now:
single-copy inverted repeats are not families, but each adds its own
length to new_bp). Tools (cd-hit-est, cd-hit-est-2d, makeblastdb,
blastn) come from PATH. Stdlib otherwise."""

import argparse
import os
import re
import subprocess
import sys

from fasta_utils import blast_subject_resolver, iter_fasta, merge_half_open, seq_id, write_fasta
from summarize_rm import read_assembly_stats

TYPES = ("tir", "helitron")
# EDTA superfamily codes -> the RepeatMasker / Dfam names in our libraries
SUPERFAMILY = {"DTA": "hAT", "DTC": "CMC", "DTH": "PIF-Harbinger", "DTM": "MULE", "DTT": "TcMar",
               "Helitron": "Helitron"}
SF_CODE = re.compile(r"(DT[ACHMT]|Helitron)")
BLAST_FMT = "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"


# ---------------------------------------------------------------- gather


def gather(args):
    os.makedirs(args.outdir, exist_ok=True)
    seen = set()
    n_in = {t: 0 for t in TYPES}
    n_out = {t: 0 for t in TYPES}
    fas = {t: open(os.path.join(args.outdir, f"{t}.candidates.fa"), "w") for t in TYPES}
    with open(os.path.join(args.outdir, "candidates.tsv"), "w") as out:
        header = None
        for tsv in args.tsv:
            fa = tsv[: -len(".tsv")] + ".fa"
            seqs = {seq_id(h).split("#", 1)[0]: (seq_id(h), s) for h, s in iter_fasta(fa)}
            with open(tsv) as fh:
                h = fh.readline()
                if header is None:
                    header = h
                    out.write(h)
                cols = h.rstrip("\n").split("\t")
                for line in fh:
                    r = dict(zip(cols, line.rstrip("\n").split("\t")))
                    t = r["type"]
                    n_in[t] += 1
                    key = (t, r["contig"], r["start"], r["end"])
                    if r["contig"] != "NA" and key in seen:
                        continue  # same element from two overlapping windows
                    seen.add(key)
                    name, seq = seqs[r["candidate"]]
                    write_fasta(fas[t], f"{args.sample}_{name}", seq)
                    out.write(line.replace(r["candidate"], f"{args.sample}_{r['candidate']}", 1))
                    n_out[t] += 1
    for fh in fas.values():
        fh.close()
    for t in TYPES:
        print(f"[dna_te] {args.sample} {t}: {n_in[t]} candidates, {n_out[t]} after dropping window-overlap duplicates")

    skipped = {}
    for path in args.timeouts:
        with open(path) as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip("\n").split("\t")
                tool, contig, s, e = f[:4]
                reason = f[4] if len(f) > 4 else "timeout"
                skipped.setdefault((tool, contig, reason), []).append((int(s), int(e)))
    with open(os.path.join(args.outdir, "skipped_windows.tsv"), "w") as out:
        out.write("tool\tcontig\tstart\tend\tbp\treason\n")
        for (tool, contig, reason), ivs in sorted(skipped.items()):
            for s, e in merge_half_open(ivs):
                out.write(f"{tool}\t{contig}\t{s}\t{e}\t{e - s}\t{reason}\n")
                print(f"[dna_te] skipped ({tool}, {reason}): {contig}:{s}-{e}")


# ---------------------------------------------------------------- pilot


def run(cmd):
    print("[dna_te] " + " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def read_clstr(path):
    """[[(name, is_rep)]] from a cd-hit .clstr; names without '>' and '...'."""
    clusters = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">Cluster"):
                clusters.append([])
                continue
            m = re.search(r">(.+?)\.\.\.", line)
            if m:
                clusters[-1].append((m.group(1), line.rstrip().endswith("*")))
    return clusters


def label_of(name):
    return name.split("#", 1)[1] if "#" in name else "NA"


def lib_names(path):
    return {seq_id(h) for h, _s in iter_fasta(path)}


def cdhit_args(args):
    return ["-c", str(args.identity), "-aS", str(args.coverage_short), "-n", str(args.word_size),
            "-G", "0", "-d", "0", "-r", "1", "-M", "0", "-T", str(args.threads)]


def match_library(args, reps_fa, lib_fa, tag, wd):
    """{rep: library family it clusters with} via cd-hit-est-2d."""
    out = os.path.join(wd, f"vs_{tag}")
    run(["cd-hit-est-2d", "-i", lib_fa, "-i2", reps_fa, "-o", out] + cdhit_args(args))
    lib = lib_names(lib_fa)
    matched = {}
    for cl in read_clstr(out + ".clstr"):
        fams = [n for n, _r in cl if n in lib]
        if fams:
            for n, _r in cl:
                if n not in lib:
                    matched[n] = fams[0]
    return matched


def read_masked(path):
    """{contig: merged [(s, e)]} 0-based, of every shared-arm .out hit."""
    ivs = {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 11 or not f[0].isdigit():
                continue
            try:
                b, e = sorted((int(f[5]), int(f[6])))
            except ValueError:
                continue
            ivs.setdefault(f[4], []).append((b - 1, e))
    return {c: merge_half_open(v) for c, v in ivs.items()}


def uncovered_bp(cover, masked):
    """bp of cover outside masked (both {contig: sorted disjoint})."""
    total = 0
    for c, ivs in cover.items():
        m = masked.get(c, [])
        j = 0
        for s, e in ivs:
            bp = e - s
            while j < len(m) and m[j][1] <= s:
                j += 1
            k = j
            while k < len(m) and m[k][0] < e:
                bp -= max(0, min(e, m[k][1]) - max(s, m[k][0]))
                k += 1
            total += bp
    return total


def superfamily(label):
    m = SF_CODE.search(label)
    return m.group(1) if m else None


def type_of(name):
    """tir / helitron from a candidate name ({sample}_{group}_{type}_{n}#label)."""
    m = re.search(r"_(tir|helitron)_\d+(#|$)", name)
    return m.group(1) if m else "NA"


# ---------------------------------------------------------------- cluster


def cluster(args):
    """cd-hit the candidates, match representatives to the libraries, and
    split the unmatched representatives into --chunks FASTA files."""
    wd = args.workdir
    os.makedirs(wd, exist_ok=True)
    skipped_bp = {}
    with open(args.skipped) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            skipped_bp[f[0]] = skipped_bp.get(f[0], 0) + int(f[4])
    shared_classes = {label_of(n) for n in lib_names(args.shared_library)}

    rows, unmatched_seqs = [], []
    with open(args.clusters_out, "w") as cl_out:
        cl_out.write("sample\ttype\tcluster_rep\tlabel\tn_members\town_match\town_class_family\t"
                     "shared_match\tshared_class_family\n")
        for t, cand in zip(TYPES, (args.tir, args.helitron)):
            n_cand = sum(1 for _ in iter_fasta(cand))
            row = {"sample": args.sample, "type": t, "candidates": n_cand, "skipped_bp": skipped_bp.get(t, 0),
                   "clusters": 0, "matched_own": 0, "matched_shared": 0, "matched": 0, "unmatched": 0,
                   "matched_classes": ".", "absent_superfamilies": "."}
            if not n_cand:
                rows.append(row)
                continue
            nr = os.path.join(wd, f"{t}.nr.fa")
            run(["cd-hit-est", "-i", cand, "-o", nr] + cdhit_args(args))
            reps = {n: len(cl) for cl in read_clstr(nr + ".clstr") for n, r in cl if r}
            own = match_library(args, nr, args.own_families, f"{t}_own", wd)
            shared = match_library(args, nr, args.shared_library, f"{t}_shared", wd)
            unmatched = [r for r in reps if r not in own and r not in shared]
            classes = {}
            for r in reps:
                fam = shared.get(r) or own.get(r)
                if fam:
                    classes[label_of(fam)] = classes.get(label_of(fam), 0) + 1
                cl_out.write(f"{args.sample}\t{t}\t{r}\t{label_of(r)}\t{reps[r]}\t{own.get(r, '.')}\t"
                             f"{label_of(own[r]) if r in own else '.'}\t{shared.get(r, '.')}\t"
                             f"{label_of(shared[r]) if r in shared else '.'}\n")
            absent = sorted({sf for r in unmatched for sf in [superfamily(label_of(r))]
                             if sf and not any(SUPERFAMILY[sf] in c for c in shared_classes)})
            keep = set(unmatched)
            unmatched_seqs += [(seq_id(h), sq) for h, sq in iter_fasta(nr) if seq_id(h) in keep]
            row.update(clusters=len(reps), matched_own=sum(1 for r in reps if r in own),
                       matched_shared=sum(1 for r in reps if r in shared),
                       matched=len(reps) - len(unmatched), unmatched=len(unmatched),
                       matched_classes=",".join(f"{c}:{n}" for c, n in
                                                sorted(classes.items(), key=lambda x: -x[1])) or ".",
                       absent_superfamilies=",".join(SUPERFAMILY[x] for x in absent) or ".")
            print(f"[dna_te] {args.sample} {t}: {row}")
            rows.append(row)

    with open(args.partial_out, "w") as out:
        out.write("\t".join(CLUSTER_COLS) + "\n")
        for r in rows:
            out.write("\t".join(str(r[c]) for c in CLUSTER_COLS) + "\n")

    # length-balanced chunks (longest first into the lightest chunk)
    loads = [[0, []] for _ in args.chunk_out]
    for name, sq in sorted(unmatched_seqs, key=lambda x: -len(x[1])):
        slot = min(loads, key=lambda x: x[0])
        slot[0] += len(sq)
        slot[1].append((name, sq))
    for path, (bp, seqs) in zip(args.chunk_out, loads):
        with open(path, "w") as fh:
            for name, sq in seqs:
                write_fasta(fh, name, sq)
    print(f"[dna_te] {len(unmatched_seqs)} unmatched representatives in {len(args.chunk_out)} chunks "
          f"({min(l[0] for l in loads)}-{max(l[0] for l in loads)} bp each)")


CLUSTER_COLS = ["sample", "type", "candidates", "skipped_bp", "clusters", "matched_own", "matched_shared",
                "matched", "unmatched", "matched_classes", "absent_superfamilies"]


# ---------------------------------------------------------------- blast_chunk


def blast_chunk(args):
    """blastn one chunk of unmatched representatives against the genome and
    reduce the HSPs to (1) per query: HSPs, distinct genomic loci, bp
    covered; (2) merged intervals per type for all queries and for the
    repeated ones (>= --min-copies loci besides the candidate's own).
    Large outputs are sorted on disk (sort -S), not held in memory."""
    os.makedirs(args.workdir, exist_ok=True)
    raw = os.path.join(args.workdir, "hsps.tsv")
    if any(True for _ in iter_fasta(args.query)):
        run(["blastn", "-task", args.task, "-query", args.query, "-db", args.db, "-evalue", str(args.evalue),
             "-max_target_seqs", str(args.max_targets), "-num_threads", str(args.threads),
             "-outfmt", "6 qseqid sseqid sstart send", "-out", raw])
    else:
        open(raw, "w").close()
    with open(args.fingerprint) as fh:
        resolve = blast_subject_resolver(l.split("\t", 1)[0] for l in fh if not l.startswith("#"))

    # query, contig, start, end (0-based half-open), sorted by query then position
    bed = os.path.join(args.workdir, "hsps.bed")
    with open(raw) as fh, open(bed, "w") as out:
        for line in fh:
            q, sub, s, e = line.rstrip("\n").split("\t")
            s, e = sorted((int(s), int(e)))
            out.write(f"{q}\t{resolve(sub)}\t{s - 1}\t{e}\n")
    os.remove(raw)
    env = dict(os.environ, LC_ALL="C")
    sort = ["sort", "-S", f"{args.sort_mem}G", "-T", args.workdir, "-t", "\t"]
    by_q = bed + ".byq"
    subprocess.run(sort + ["-k1,1", "-k2,2", "-k3,3n", "-o", by_q, bed], check=True, env=env)

    queries = {}  # q -> [hsps, loci, cover_bp]
    cur_q, cur_c, ivs, n_hsp, loci = None, None, [], 0, []
    with open(by_q) as fh:
        for line in fh:
            q, c, s, e = line.rstrip("\n").split("\t")
            if q != cur_q:
                if cur_q is not None:
                    loci += [(cur_c, iv) for iv in merge_half_open(ivs)]
                    queries[cur_q] = [n_hsp, len(loci), sum(iv[1] - iv[0] for _c, iv in loci)]
                cur_q, cur_c, ivs, n_hsp, loci = q, c, [], 0, []
            elif c != cur_c:
                loci += [(cur_c, iv) for iv in merge_half_open(ivs)]
                cur_c, ivs = c, []
            ivs.append((int(s), int(e)))
            n_hsp += 1
        if cur_q is not None:
            loci += [(cur_c, iv) for iv in merge_half_open(ivs)]
            queries[cur_q] = [n_hsp, len(loci), sum(iv[1] - iv[0] for _c, iv in loci)]
    os.remove(by_q)
    for name, _sq in iter_fasta(args.query):
        queries.setdefault(seq_id(name), [0, 0, 0])
    repeated = {q for q, v in queries.items() if v[1] >= args.min_copies + 1}

    with open(args.queries_out, "w") as out:
        out.write("query\ttype\thsps\tloci\tcover_bp\trepeated\n")
        for q, (n, l, bp) in sorted(queries.items()):
            out.write(f"{q}\t{type_of(q)}\t{n}\t{l}\t{bp}\t{q in repeated}\n")

    # merged intervals per (type, set), sorted by contig on disk
    tagged = bed + ".tag"
    with open(bed) as fh, open(tagged, "w") as out:
        for line in fh:
            q, c, s, e = line.rstrip("\n").split("\t")
            t = type_of(q)
            out.write(f"{t}\tall\t{c}\t{s}\t{e}\n")
            if q in repeated:
                out.write(f"{t}\trepeated\t{c}\t{s}\t{e}\n")
    os.remove(bed)
    by_c = tagged + ".sorted"
    subprocess.run(sort + ["-k1,1", "-k2,2", "-k3,3", "-k4,4n", "-o", by_c, tagged], check=True, env=env)
    os.remove(tagged)
    with open(by_c) as fh, open(args.bed_out, "w") as out:
        key, cs, ce = None, None, None
        for line in fh:
            t, st, c, s, e = line.rstrip("\n").split("\t")
            s, e = int(s), int(e)
            if key == (t, st, c) and s <= ce:
                ce = max(ce, e)
                continue
            if key is not None:
                out.write(f"{key[0]}\t{key[1]}\t{key[2]}\t{cs}\t{ce}\n")
            key, cs, ce = (t, st, c), s, e
        if key is not None:
            out.write(f"{key[0]}\t{key[1]}\t{key[2]}\t{cs}\t{ce}\n")
    os.remove(by_c)
    print(f"[dna_te] {args.query}: {len(queries)} queries, {len(repeated)} repeated "
          f"(>= {args.min_copies} copies besides their own locus)")


# ---------------------------------------------------------------- finalize


def finalize(args):
    """Chunk summaries -> pilot.tsv (per type) and pilot_clusters.tsv with copies."""
    _total, non_n = read_assembly_stats(args.assembly_stats)
    masked = read_masked(args.out_file)
    ivs = {}
    for path in args.beds:
        with open(path) as fh:
            for line in fh:
                t, st, c, s, e = line.rstrip("\n").split("\t")
                ivs.setdefault((t, st), {}).setdefault(c, []).append((int(s), int(e)))
    cover = {k: {c: merge_half_open(v) for c, v in d.items()} for k, d in ivs.items()}
    qinfo = {}
    for path in args.queries:
        with open(path) as fh:
            fh.readline()
            for line in fh:
                q, t, n, l, bp, rep = line.rstrip("\n").split("\t")
                qinfo[q] = (n, l, bp, rep)

    def bp_of(t, st):
        d = cover.get((t, st), {})
        return sum(e - s for v in d.values() for s, e in v), uncovered_bp(d, masked)

    with open(args.partial) as fh, open(args.out, "w") as out:
        cols = fh.readline().rstrip("\n").split("\t")
        out.write("\t".join(PILOT_COLS) + "\n")
        for line in fh:
            r = dict(zip(cols, line.rstrip("\n").split("\t")))
            t = r["type"]
            r["unmatched_cover_bp"], r["unmatched_new_bp"] = bp_of(t, "all")
            _rc, r["new_bp_repeated"] = bp_of(t, "repeated")
            r["unmatched_repeated"] = sum(1 for q, v in qinfo.items() if type_of(q) == t and v[3] == "True")
            for k, pct in (("unmatched_new_bp", "new_pct_non_n"), ("new_bp_repeated", "new_pct_repeated")):
                r[pct] = f"{100.0 * r[k] / non_n:.4f}" if non_n else "NA"
            out.write("\t".join(str(r[c]) for c in PILOT_COLS) + "\n")
            print(f"[dna_te] {r['sample']} {t}: {r}")

    with open(args.clusters_in) as fh, open(args.clusters_out, "w") as out:
        out.write(fh.readline().rstrip("\n") + "\tgenome_hsps\tgenome_loci\tgenome_cover_bp\trepeated\n")
        for line in fh:
            rep = line.split("\t")[2]
            n, l, bp, rp = qinfo.get(rep, (".", ".", ".", "."))
            out.write(line.rstrip("\n") + f"\t{n}\t{l}\t{bp}\t{rp}\n")


PILOT_COLS = ["sample", "type", "candidates", "skipped_bp", "clusters", "matched_own", "matched_shared", "matched",
              "unmatched", "unmatched_cover_bp", "unmatched_new_bp", "new_pct_non_n", "unmatched_repeated",
              "new_bp_repeated", "new_pct_repeated", "matched_classes", "absent_superfamilies"]


# ---------------------------------------------------------------- combine


def combine(args):
    with open(args.out, "w") as out:
        out.write(f"# decision rule: integrate (Stage B2) if the unmatched clusters add "
                  f">= {args.min_new_pct_non_n} % of non-N bp the shared library does not mask yet "
                  f"(new_pct_non_n), or carry a superfamily absent from the shared library "
                  f"(absent_superfamilies). new bp is a blastn lower bound.\n")
        out.write("\t".join(PILOT_COLS + ["decision"]) + "\n")
        for path in args.chunks:
            with open(path) as fh:
                cols = fh.readline().rstrip("\n").split("\t")
                for line in fh:
                    r = dict(zip(cols, line.rstrip("\n").split("\t")))
                    pct = float(r["new_pct_non_n"]) if r["new_pct_non_n"] != "NA" else 0.0
                    why = []
                    if pct >= args.min_new_pct_non_n:
                        why.append(f"adds {pct:.2f}%")
                    if r["absent_superfamilies"] != ".":
                        why.append(f"new superfamily {r['absent_superfamilies']}")
                    dec = "integrate: " + "; ".join(why) if why else "skip"
                    out.write("\t".join(r[c] for c in PILOT_COLS) + f"\t{dec}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    g = sub.add_parser("gather")
    g.add_argument("--sample", required=True)
    g.add_argument("--tsv", nargs="+", required=True, help="dna_te_windows.py --out-tsv files (.fa next to each)")
    g.add_argument("--timeouts", nargs="*", default=[])
    g.add_argument("--outdir", required=True)

    p = sub.add_parser("cluster", help="cd-hit + library matching; unmatched reps split into chunks")
    p.add_argument("--sample", required=True)
    p.add_argument("--tir", required=True)
    p.add_argument("--helitron", required=True)
    p.add_argument("--own-families", required=True)
    p.add_argument("--shared-library", required=True)
    p.add_argument("--skipped", required=True, help="skipped_windows.tsv from gather (bp not scanned per type)")
    p.add_argument("--identity", type=float, default=0.8)
    p.add_argument("--coverage-short", type=float, default=0.8)
    p.add_argument("--word-size", type=int, default=5)
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--workdir", required=True)
    p.add_argument("--partial-out", required=True)
    p.add_argument("--clusters-out", required=True)
    p.add_argument("--chunk-out", nargs="+", required=True, help="one FASTA per blast chunk")

    b = sub.add_parser("blast_chunk", help="blastn one chunk against the genome; compact summaries")
    b.add_argument("--query", required=True)
    b.add_argument("--db", required=True)
    b.add_argument("--fingerprint", required=True, help="genome fingerprint.tsv (sequence names)")
    b.add_argument("--task", default="dc-megablast")
    b.add_argument("--evalue", type=float, default=1e-10)
    b.add_argument("--max-targets", type=int, required=True)
    b.add_argument("--min-copies", type=int, default=3, help="loci besides the candidate's own for 'repeated'")
    b.add_argument("--threads", type=int, default=1)
    b.add_argument("--sort-mem", type=int, default=2, help="GB for sort -S")
    b.add_argument("--workdir", required=True)
    b.add_argument("--queries-out", required=True)
    b.add_argument("--bed-out", required=True)

    f = sub.add_parser("finalize", help="chunk summaries -> pilot.tsv and pilot_clusters.tsv")
    f.add_argument("--partial", required=True)
    f.add_argument("--clusters-in", required=True)
    f.add_argument("--queries", nargs="+", required=True)
    f.add_argument("--beds", nargs="+", required=True)
    f.add_argument("--out-file", required=True, help="shared-arm .out (bp already masked)")
    f.add_argument("--assembly-stats", required=True)
    f.add_argument("--out", required=True)
    f.add_argument("--clusters-out", required=True)

    c = sub.add_parser("combine")
    c.add_argument("--chunks", nargs="+", required=True)
    c.add_argument("--min-new-pct-non-n", type=float, default=0.5)
    c.add_argument("--out", required=True)

    args = ap.parse_args()
    {"gather": gather, "cluster": cluster, "blast_chunk": blast_chunk, "finalize": finalize,
     "combine": combine}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
