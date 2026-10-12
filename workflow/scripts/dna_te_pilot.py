#!/usr/bin/env python3
"""Structural TIR / Helitron pilot (README "Structural TIR / Helitron
pilot"). Report-only: the library is not changed.

  gather   per-group candidates from dna_te_windows.py -> one FASTA per type
           (candidates reported twice by overlapping windows dropped),
           candidates.tsv, and skipped_windows.tsv (timed-out or failed
           windows, merged per tool and reason, as ltr/skipped_windows.tsv).
  cluster  per sample and type: cd-hit-est the candidates (library.cdhit
           settings, both strands); cd-hit-est-2d the cluster
           representatives against the sample's own -families.fa and the
           shared library. The unmatched ones are blastn'd against the shared
           library and every segment matching a family (>= --mask-min-id %)
           is hard-masked: a representative with < --novel-min-bp or
           < --novel-min-frac left is a *fragment* of its largest library
           match (e.g. a CR1 piece TIR-Learner called DTC); the rest are
           *novel* and go, masked, into length-balanced chunks. Candidates
           covering >= --unknown-min-cov of an Unknown family's consensus
           are listed in unknown_overlaps.tsv (structural evidence for it).
  blast_chunk  blastn (dc-megablast) one chunk against the genome; reduce
           the HSPs on disk to per-query copies and merged intervals.
  finalize per sample: chunk summaries -> one row per type (pilot.tsv), the
           per-cluster table with genome copies, novel_repeated.fa.
  targets  one TEsorter FASTA: repeated novel representatives, Unknown
           families the candidates overlap, top Unknown families of
           curation_candidates.tsv.
  report   pilot rows + TEsorter -> summary/dna_te_pilot.tsv (decision rule
           on its first line) and dna_te_unknown_proposals.tsv
           (classify.curated_families format, review before use).

The bp the novel clusters would mask is a blastn proxy: the union of their
dc-megablast HSPs on the genome, less sensitive than RepeatMasker -s (old,
diverged copies are missed), so a lower bound. `novel_new_bp` is the part
the shared-arm .out does not mask already; `*_repeated` counts only novel
representatives with >= --min-copies genomic loci besides their own, and
the decision uses new_pct_repeated. Tools (cd-hit-est, cd-hit-est-2d, makeblastdb,
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

# TIR-Learner / HelitronScanner superfamily -> RepeatMasker Class/Family
EDTA_TO_RM = {"DTA": "DNA/hAT", "DTC": "DNA/CMC-EnSpm", "DTH": "DNA/PIF-Harbinger", "DTM": "DNA/MULE-MuDR",
              "DTT": "DNA/TcMar", "Helitron": "RC/Helitron"}


def library_labels(path, curated_paths):
    """{bare family name: Class/Family} from the library headers, overridden
    by the curated / applied cross-check tables (first table listing a family
    wins, as in summarize)."""
    labels = {seq_id(h).split("#", 1)[0]: label_of(seq_id(h)) for h, _s in iter_fasta(path)}
    if curated_paths:
        from family_groups import curated_classes
        seen = set()
        for cp in curated_paths:
            for fam, cf in curated_classes(cp).items():
                if fam not in seen:
                    labels[fam] = cf
                    seen.add(fam)
    return labels


def mask_against_library(args, seqs, wd):
    """seqs: {name: seq}. blastn them against the shared library; return
    {name: (masked_ivs, per_family)} with masked_ivs the merged 0-based query
    intervals covered by library HSPs at >= --mask-min-id, and per_family
    {family: [query_bp, family_cov, best_pident]}."""
    fa = os.path.join(wd, "unmatched.fa")
    with open(fa, "w") as fh:
        for n, sq in seqs.items():
            write_fasta(fh, n, sq)
    db = os.path.join(wd, "libdb")
    run(["makeblastdb", "-in", args.shared_library, "-dbtype", "nucl", "-out", db])
    out = os.path.join(wd, "unmatched_vs_library.tsv")
    mt = ["-mt_mode", "1"] if args.threads > 1 else []
    run(["blastn", "-task", "dc-megablast", "-query", fa, "-db", db, "-evalue", str(args.evalue),
         "-max_target_seqs", "1000", "-num_threads", str(args.threads)] + mt +
        ["-outfmt", "6 qseqid sseqid pident qstart qend sstart send slen", "-out", out])
    q_ivs, fam_q, fam_s, fam_id, slen = {}, {}, {}, {}, {}
    with open(out) as fh:
        for line in fh:
            q, sub, pid, qs, qe, ss, se, sl = line.rstrip("\n").split("\t")
            if float(pid) < args.mask_min_id:
                continue
            fam = sub.split("#", 1)[0]
            qs, qe = int(qs) - 1, int(qe)
            ss, se = sorted((int(ss), int(se)))
            q_ivs.setdefault(q, []).append((qs, qe))
            fam_q.setdefault((q, fam), []).append((qs, qe))
            fam_s.setdefault((q, fam), []).append((ss - 1, se))
            fam_id[(q, fam)] = max(fam_id.get((q, fam), 0.0), float(pid))
            slen[fam] = int(sl)
    res = {n: ([], {}) for n in seqs}
    for q, ivs in q_ivs.items():
        res[q] = (merge_half_open(ivs), res[q][1])
    for (q, fam), ivs in fam_q.items():
        qbp = sum(e - s for s, e in merge_half_open(ivs))
        sbp = sum(e - s for s, e in merge_half_open(fam_s[(q, fam)]))
        res[q][1][fam] = [qbp, sbp / slen[fam] if slen[fam] else 0.0, fam_id[(q, fam)]]
    return res


def hard_mask(seq, ivs):
    s = list(seq)
    for a, b in ivs:
        s[a:b] = "N" * (b - a)
    return "".join(s)


def cluster(args):
    """cd-hit the candidates, match representatives to the libraries, mask the
    unmatched ones against the shared library, and split the novel ones
    (masked) into --chunks FASTA files."""
    wd = args.workdir
    os.makedirs(wd, exist_ok=True)
    skipped_bp = {}
    with open(args.skipped) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            skipped_bp[f[0]] = skipped_bp.get(f[0], 0) + int(f[4])
    # header labels only: the current labels (curated / cross-check tables) are
    # applied later, in targets / report, so relabelling never reruns this step
    labels = library_labels(args.shared_library, [])
    shared_classes = set(labels.values())

    rows, table, unmatched_seqs = {}, [], {}
    for t, cand in zip(TYPES, (args.tir, args.helitron)):
        n_cand = sum(1 for _ in iter_fasta(cand))
        row = {"sample": args.sample, "type": t, "candidates": n_cand, "skipped_bp": skipped_bp.get(t, 0),
               "clusters": 0, "matched_own": 0, "matched_shared": 0, "matched": 0, "unmatched": 0,
               "fragments": 0, "novel": 0, "matched_classes": ".", "fragment_classes": ".",
               "absent_superfamilies": "."}
        rows[t] = row
        if not n_cand:
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
            table.append({"t": t, "rep": r, "n": reps[r], "own": own.get(r), "shared": shared.get(r)})
        keep = set(unmatched)
        unmatched_seqs.update({seq_id(h): sq for h, sq in iter_fasta(nr) if seq_id(h) in keep})
        row.update(clusters=len(reps), matched_own=sum(1 for r in reps if r in own),
                   matched_shared=sum(1 for r in reps if r in shared),
                   matched=len(reps) - len(unmatched), unmatched=len(unmatched),
                   matched_classes=",".join(f"{c}:{n}" for c, n in
                                            sorted(classes.items(), key=lambda x: -x[1])) or ".")

    masks = mask_against_library(args, unmatched_seqs, wd) if unmatched_seqs else {}
    status, novel_seqs, frag_classes = {}, [], {t: {} for t in TYPES}
    with open(args.overlaps_out, "w") as ov:
        # every family a candidate covers >= unknown_min_cov of (header label);
        # targets / report keep the ones whose *current* label is Unknown
        ov.write("sample\ttype\tcandidate\tcandidate_label\tfamily\tfamily_label\tfamily_cov\tquery_bp\tpident\n")
        for name, sq in unmatched_seqs.items():
            t = type_of(name)
            ivs, per_fam = masks.get(name, ([], {}))
            masked = sum(e - s for s, e in ivs)
            novel_bp = len(sq) - masked
            top = max(per_fam.items(), key=lambda x: x[1][0])[0] if per_fam else None
            if novel_bp < args.novel_min_bp or novel_bp < args.novel_min_frac * len(sq):
                lab = labels.get(top, "NA")
                status[name] = ("fragment", masked / len(sq), f"{top}#{lab}")
                frag_classes[t][lab] = frag_classes[t].get(lab, 0) + 1
                rows[t]["fragments"] += 1
            else:
                status[name] = ("novel", masked / len(sq), f"{top}#{labels.get(top, 'NA')}" if top else ".")
                novel_seqs.append((name, hard_mask(sq, ivs)))
                rows[t]["novel"] += 1
            for fam, (qbp, fcov, pid) in per_fam.items():
                if fcov >= args.unknown_min_cov:
                    ov.write(f"{args.sample}\t{t}\t{name}\t{label_of(name)}\t{fam}\t{labels[fam]}\t"
                             f"{fcov:.3f}\t{qbp}\t{pid:.1f}\n")

    for t, row in rows.items():
        novel_names = [n for n, (st, _f, _x) in status.items() if st == "novel" and type_of(n) == t]
        absent = sorted({sf for n in novel_names for sf in [superfamily(label_of(n))]
                         if sf and not any(SUPERFAMILY[sf] in c for c in shared_classes)})
        row["absent_superfamilies"] = ",".join(SUPERFAMILY[x] for x in absent) or "."
        row["fragment_classes"] = ",".join(f"{c}:{n}" for c, n in sorted(frag_classes[t].items(),
                                                                         key=lambda x: -x[1])) or "."
        print(f"[dna_te] {args.sample} {t}: {row}")

    with open(args.partial_out, "w") as out:
        out.write("\t".join(CLUSTER_COLS) + "\n")
        for t in TYPES:
            out.write("\t".join(str(rows[t][c]) for c in CLUSTER_COLS) + "\n")

    with open(args.clusters_out, "w") as cl_out:
        cl_out.write("sample\ttype\tcluster_rep\tlabel\tn_members\town_match\town_class_family\t"
                     "shared_match\tshared_class_family\tstatus\tlibrary_masked_frac\tlargest_library_overlap\n")
        for r in table:
            own, sh = r["own"], r["shared"]
            st, mf, top = status.get(r["rep"], ("matched", None, "."))
            cl_out.write(f"{args.sample}\t{r['t']}\t{r['rep']}\t{label_of(r['rep'])}\t{r['n']}\t{own or '.'}\t"
                         f"{label_of(own) if own else '.'}\t{sh or '.'}\t{label_of(sh) if sh else '.'}\t"
                         f"{st}\t{'.' if mf is None else f'{mf:.3f}'}\t{top}\n")

    # length-balanced chunks of the novel, masked representatives
    loads = [[0, []] for _ in args.chunk_out]
    for name, sq in sorted(novel_seqs, key=lambda x: -len(x[1])):
        slot = min(loads, key=lambda x: x[0])
        slot[0] += len(sq)
        slot[1].append((name, sq))
    for path, (bp, seqs) in zip(args.chunk_out, loads):
        with open(path, "w") as fh:
            for name, sq in seqs:
                write_fasta(fh, name, sq)
    print(f"[dna_te] {len(unmatched_seqs)} unmatched representatives: "
          f"{sum(1 for v in status.values() if v[0] == 'fragment')} fragments of library families, "
          f"{len(novel_seqs)} novel (masked) in {len(args.chunk_out)} chunks")


CLUSTER_COLS = ["sample", "type", "candidates", "skipped_bp", "clusters", "matched_own", "matched_shared",
                "matched", "unmatched", "fragments", "novel", "matched_classes", "fragment_classes",
                "absent_superfamilies"]


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
        # -mt_mode 1: threads split the *queries*. The default (0) splits the
        # database, which only parallelises the seed scan; the gapped
        # traceback of a repeat fragment's tens of thousands of hits then
        # runs on one thread while the others idle.
        mt = ["-mt_mode", "1"] if args.threads > 1 else []
        run(["blastn", "-task", args.task, "-query", args.query, "-db", args.db, "-evalue", str(args.evalue),
             "-max_target_seqs", str(args.max_targets), "-num_threads", str(args.threads)] + mt +
            ["-outfmt", "6 qseqid sseqid sstart send", "-out", raw])
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
    """Chunk summaries -> pilot.tsv (per type), pilot_clusters.tsv with genome
    copies, and novel_repeated.fa (the repeated novel representatives, masked,
    for TEsorter)."""
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
            r["novel_cover_bp"], r["novel_new_bp"] = bp_of(t, "all")
            _rc, r["new_bp_repeated"] = bp_of(t, "repeated")
            r["novel_repeated"] = sum(1 for q, v in qinfo.items() if type_of(q) == t and v[3] == "True")
            for k, pct in (("novel_new_bp", "new_pct_novel"), ("new_bp_repeated", "new_pct_repeated")):
                r[pct] = f"{100.0 * r[k] / non_n:.4f}" if non_n else "NA"
            out.write("\t".join(str(r[c]) for c in PILOT_COLS) + "\n")
            print(f"[dna_te] {r['sample']} {t}: {r}")

    with open(args.clusters_in) as fh, open(args.clusters_out, "w") as out:
        out.write(fh.readline().rstrip("\n") + "\tgenome_hsps\tgenome_loci\tgenome_cover_bp\trepeated\n")
        for line in fh:
            rep = line.split("\t")[2]
            n, l, bp, rp = qinfo.get(rep, (".", ".", ".", "."))
            out.write(line.rstrip("\n") + f"\t{n}\t{l}\t{bp}\t{rp}\n")

    n = 0
    with open(args.novel_out, "w") as out:
        for path in args.chunks:
            for h, sq in iter_fasta(path):
                if qinfo.get(seq_id(h), (0, 0, 0, "False"))[3] == "True":
                    write_fasta(out, seq_id(h), sq)
                    n += 1
    print(f"[dna_te] {n} repeated novel representatives -> {args.novel_out}")


PILOT_COLS = ["sample", "type", "candidates", "skipped_bp", "clusters", "matched_own", "matched_shared", "matched",
              "unmatched", "fragments", "novel", "novel_cover_bp", "novel_new_bp", "new_pct_novel",
              "novel_repeated", "new_bp_repeated", "new_pct_repeated", "matched_classes", "fragment_classes",
              "absent_superfamilies"]


# ---------------------------------------------------------------- targets


def targets(args):
    """One FASTA for TEsorter: repeated novel representatives (all samples),
    Unknown families TIR/Helitron candidates overlap, and the top Unknown
    families of curation_candidates.tsv (flag `unknown`, not yet curated)."""
    src = {}  # bare name -> {"sources": set, "samples": set, "type": str}

    def add(name, source, sample, t="."):
        e = src.setdefault(name, {"sources": set(), "samples": set(), "type": t})
        e["sources"].add(source)
        e["samples"].add(sample)
    seqs = {}
    for path in args.novel:
        for h, sq in iter_fasta(path):
            name = seq_id(h).split("#", 1)[0]
            seqs[name] = sq
            add(name, "novel_repeated", name.split("_", 1)[0], type_of(seq_id(h)))
    labels = library_labels(args.library, args.curated_table)
    for path in args.overlaps:
        with open(path) as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if labels.get(f[4], "").split("/")[0] == "Unknown":
                    add(f[4], "unknown_overlap", f[0])
    with open(args.curation) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            if "unknown" in r.get("flags", "").split(",") and r.get("curated", ".") == ".":
                add(r["family"], "top_unknown", "curation_candidates")
    lib = {seq_id(h).split("#", 1)[0]: sq for h, sq in iter_fasta(args.library)}
    missing = 0
    with open(args.fa_out, "w") as fa, open(args.tsv_out, "w") as tsv:
        tsv.write("name\ttype\tsources\tsamples\tlength\n")
        for name, e in sorted(src.items()):
            sq = seqs.get(name) or lib.get(name)
            if sq is None:
                missing += 1
                continue
            write_fasta(fa, name, sq)
            tsv.write(f"{name}\t{e['type']}\t{','.join(sorted(e['sources']))}\t"
                      f"{','.join(sorted(e['samples']))}\t{len(sq)}\n")
    by = {}
    for e in src.values():
        for x in e["sources"]:
            by[x] = by.get(x, 0) + 1
    print(f"[dna_te] TEsorter targets: {len(src) - missing} ({by}); {missing} names not found in the library")


# ---------------------------------------------------------------- report


def report(args):
    """Per-sample pilot rows + TEsorter calls -> summary/dna_te_pilot.tsv and
    the Unknown-family label proposals (classify.curated_families format)."""
    from reclassify_unknown import domain_evidence, read_tesorter_cls
    rex, gydb = read_tesorter_cls(args.rexdb), read_tesorter_cls(args.gydb)
    tg = {}
    with open(args.targets) as fh:
        fh.readline()
        for line in fh:
            name, t, sources, samples, _l = line.rstrip("\n").split("\t")
            tg[name] = (t, set(sources.split(",")), set(samples.split(",")))
    dom = {n: domain_evidence(n, rex, gydb, args.min_domains) for n in tg}
    lib_classes = {label_of(seq_id(h)) for h, _s in iter_fasta(args.library)}

    with open(args.out, "w") as out:
        out.write(f"# decision rule: integrate (Stage B2) if the repeated novel clusters (>= min_copies genomic "
                  f"copies; segments matching a library family masked out first) add >= {args.min_new_pct_non_n} % "
                  f"of non-N bp the shared library does not mask yet (new_pct_repeated), or carry a TE domain of a "
                  f"superfamily absent from the shared library. `review` = enough bp but no repeated novel "
                  f"cluster carries a TE domain. new bp is a blastn lower bound.\n")
        cols_out = None
        for path in args.chunks:
            with open(path) as fh:
                cols = fh.readline().rstrip("\n").split("\t")
                if cols_out is None:
                    cols_out = cols + ["repeated_with_domain", "domain_classes", "decision"]
                    out.write("\t".join(cols_out) + "\n")
                for line in fh:
                    r = dict(zip(cols, line.rstrip("\n").split("\t")))
                    mine = [n for n, (t, srcs, smp) in tg.items() if "novel_repeated" in srcs
                            and t == r["type"] and r["sample"] in smp]
                    calls = [dom[n][0] for n in mine if dom[n][0]]
                    classes = {}
                    for c in calls:
                        classes[c] = classes.get(c, 0) + 1
                    te_new = sorted({c for c in calls if c.split("/")[0] in ("DNA", "RC") and "/" in c
                                     and not any(c == x or x.startswith(c + "-") or c.startswith(x + "-")
                                                 for x in lib_classes)})
                    pct = float(r["new_pct_repeated"]) if r["new_pct_repeated"] != "NA" else 0.0
                    why = []
                    if pct >= args.min_new_pct_non_n:
                        why.append(f"adds {pct:.2f}% ({len(calls)} of {len(mine)} repeated novel clusters "
                                   f"with a TE domain)")
                    if te_new:
                        why.append(f"domain-supported superfamily absent from the library: {','.join(te_new)}")
                    if not why:
                        dec = "skip"
                    elif calls:
                        dec = "integrate: " + "; ".join(why)
                    else:
                        dec = "review: " + "; ".join(why)
                    r.update(repeated_with_domain=len(calls),
                             domain_classes=",".join(f"{c}:{n}" for c, n in
                                                     sorted(classes.items(), key=lambda x: -x[1])) or ".",
                             decision=dec)
                    out.write("\t".join(str(r[c]) for c in cols_out) + "\n")

    # Unknown-family proposals: TEsorter domains and/or TIR/Helitron structure
    struct = {}
    for path in args.overlaps:
        with open(path) as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip("\n").split("\t")
                sf = superfamily(f[3])
                struct.setdefault(f[4], []).append((f[0], f[2], EDTA_TO_RM.get(sf, "NA"), float(f[6])))
    n_rows = {}
    with open(args.proposals_out, "w") as out:
        out.write("group\tfamily\tclass_family\tpart\tevidence\tnote\tconfidence\n")
        for fam, (t, srcs, smp) in sorted(tg.items()):
            if "novel_repeated" in srcs and len(srcs) == 1:
                continue
            d, conflict = dom.get(fam, (None, ""))
            st = struct.get(fam, [])
            s_classes = sorted({x[2] for x in st if x[2] != "NA"})
            s_note = (f"{len(st)} TIR/Helitron candidates cover >= {args.unknown_min_cov} of it "
                      f"({','.join(s_classes) or 'NA'}; best cov {max(x[3] for x in st):.2f})") if st else ""
            dom_note = f"TEsorter {d}" if d else (f"TEsorter conflict: {conflict}" if conflict else "")
            if d and st:
                agree = any(d.split("/")[0] == c.split("/")[0] for c in s_classes)
                if agree:
                    cf, conf, ev = d, "high", "tesorter+structure"
                else:
                    cf, conf, ev = ".", "conflict", "tesorter_vs_structure"
            elif d:
                cf, conf, ev = d, "medium", "tesorter"
            elif st:
                root = sorted({c.split("/")[0] for c in s_classes}) or ["DNA"]
                cf = "RC/Helitron" if root == ["RC"] else ("DNA" if root == ["DNA"] else ".")
                conf, ev = ("low", "structure") if cf != "." else ("conflict", "structure_mixed")
            else:
                continue
            note = "; ".join(x for x in (dom_note, s_note, f"sources {','.join(sorted(srcs))}") if x)
            out.write(f"{fam}\t{fam}\t{cf}\t.\tdna_te_pilot:{ev}\t{note}\t{conf}\n")
            n_rows[conf] = n_rows.get(conf, 0) + 1
    print(f"[dna_te] Unknown-family proposals: {n_rows}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    g = sub.add_parser("gather")
    g.add_argument("--sample", required=True)
    g.add_argument("--tsv", nargs="+", required=True, help="dna_te_windows.py --out-tsv files (.fa next to each)")
    g.add_argument("--timeouts", nargs="*", default=[])
    g.add_argument("--outdir", required=True)

    p = sub.add_parser("cluster", help="cd-hit + library matching + masking; novel reps split into chunks")
    p.add_argument("--sample", required=True)
    p.add_argument("--tir", required=True)
    p.add_argument("--helitron", required=True)
    p.add_argument("--own-families", required=True)
    p.add_argument("--shared-library", required=True)
    p.add_argument("--skipped", required=True, help="skipped_windows.tsv from gather (bp not scanned per type)")
    p.add_argument("--identity", type=float, default=0.8)
    p.add_argument("--coverage-short", type=float, default=0.8)
    p.add_argument("--word-size", type=int, default=5)
    p.add_argument("--evalue", type=float, default=1e-10)
    p.add_argument("--mask-min-id", type=float, default=70.0, help="library HSPs at >= this %% identity mask")
    p.add_argument("--novel-min-bp", type=int, default=200)
    p.add_argument("--novel-min-frac", type=float, default=0.5)
    p.add_argument("--unknown-min-cov", type=float, default=0.5,
                   help="a candidate covering >= this of an Unknown family's consensus is recorded")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--workdir", required=True)
    p.add_argument("--partial-out", required=True)
    p.add_argument("--clusters-out", required=True)
    p.add_argument("--overlaps-out", required=True)
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

    f = sub.add_parser("finalize", help="chunk summaries -> pilot.tsv, pilot_clusters.tsv, novel_repeated.fa")
    f.add_argument("--partial", required=True)
    f.add_argument("--clusters-in", required=True)
    f.add_argument("--queries", nargs="+", required=True)
    f.add_argument("--beds", nargs="+", required=True)
    f.add_argument("--chunks", nargs="+", required=True, help="the chunk FASTAs (masked novel reps)")
    f.add_argument("--out-file", required=True, help="shared-arm .out (bp already masked)")
    f.add_argument("--assembly-stats", required=True)
    f.add_argument("--out", required=True)
    f.add_argument("--clusters-out", required=True)
    f.add_argument("--novel-out", required=True)

    tg = sub.add_parser("targets", help="TEsorter targets: novel repeated + overlapped / top Unknown families")
    tg.add_argument("--novel", nargs="+", required=True)
    tg.add_argument("--overlaps", nargs="+", required=True)
    tg.add_argument("--curation", required=True, help="curation_candidates.tsv")
    tg.add_argument("--library", required=True)
    tg.add_argument("--curated-table", nargs="*", default=[],
                    help="curated / applied cross-check tables: current labels (first listing wins)")
    tg.add_argument("--fa-out", required=True)
    tg.add_argument("--tsv-out", required=True)

    c = sub.add_parser("report", help="pilot rows + TEsorter -> summary table and Unknown proposals")
    c.add_argument("--chunks", nargs="+", required=True, help="per-sample pilot.tsv")
    c.add_argument("--targets", required=True)
    c.add_argument("--rexdb", required=True)
    c.add_argument("--gydb", required=True)
    c.add_argument("--overlaps", nargs="+", required=True)
    c.add_argument("--library", required=True)
    c.add_argument("--min-domains", type=int, default=1)
    c.add_argument("--unknown-min-cov", type=float, default=0.5)
    c.add_argument("--min-new-pct-non-n", type=float, default=0.5)
    c.add_argument("--out", required=True)
    c.add_argument("--proposals-out", required=True)

    args = ap.parse_args()
    {"gather": gather, "cluster": cluster, "blast_chunk": blast_chunk, "finalize": finalize,
     "targets": targets, "report": report}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
