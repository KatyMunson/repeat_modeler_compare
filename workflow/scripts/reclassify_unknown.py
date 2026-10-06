#!/usr/bin/env python3
"""Reclassify the library's Unknown families from three cheap, independent
screens of their consensus sequences:

  TEsorter   TE protein domains (REXdb metazoa, GyDB) -> order / superfamily
  Rfam       structured RNA (cmscan --cut_ga)         -> rRNA / tRNA / snRNA / ...
  DIAMOND    host proteins (TE proteins removed)      -> Other/host_gene

Subcommands:
  extract   write the library's records, header = bare family name (TEsorter
            and cmscan mangle '#'), plus a family -> RepeatModeler class table.
            --keep-prefix limits it to de novo families (the pipeline passes
            each sample's prefix, so curated Dfam entries are left alone);
            --unknown-only keeps only Unknown ones.
  merge     combine the screens into one reclassification table for the
            Unknown families
  verify    report-only check of the CLASSIFIED families: does the same
            independent evidence agree with RepeatClassifier's label?
            bp-weighted with family_tandem.tsv tables (see cmd_verify)

Precedence (merge), first match wins:
  1. TEsorter call at order level or finer, REXdb and GyDB not disagreeing
     on order, and >= --min-domains domains           evidence = domain
  2. Rfam hit covering >= --min-rfam-cov of the consensus  evidence = rfam
  3. DIAMOND host hit covering >= --min-host-cov, and the family has
     <= --host-max-copies .out hits in every sample     evidence = host_protein
  Disagreements (REXdb vs GyDB order; a domain call AND a host hit; a domain
  call AND an Rfam hit) leave the family Unknown with `conflict` recorded.
A domain call beats Rfam and host hits only when they don't also qualify:
both qualifying is a conflict, not a win, because a host gene or RNA
carrying a TE domain (a domesticated TE, a TE-derived exon) is ambiguous.

Host proteins come from genome annotations, which turn many TE open reading
frames into "uncharacterized LOC" genes that carry no recognisable TE domain
(so the keyword and TEsorter filters miss them). A family with hundreds of
genomic copies matching such a protein is a TE, not a host gene: above
--host-max-copies (needs --family-tandem for copy counts) the host hit is
ignored for labelling, the family stays Unknown with
note = host_annotated_te_orf, and verify reports the status of the same
name instead of a host_protein disagreement.

Only families that change are relabelled, but every Unknown family gets a row
and diamond_best / rfam_best report the best hit even below the coverage cut
(a tRNA-headed SINE shows a partial tRNA hit), so the table also documents
what each screen saw. summarize_rm.py applies
the table (--reclass-table) to .out/.align class labels; nothing is remasked.
Stdlib only."""

import argparse
import os
import re
import sys

from fasta_utils import iter_fasta, write_fasta

# TEsorter order -> RepeatMasker class. Orders not listed fall back to
# leaving the family Unknown (recorded in the table, not relabelled).
ORDER_TO_RM = {
    "LTR": "LTR",
    "LINE": "LINE",
    "SINE": "SINE",
    "DIRS": "LTR/DIRS",
    "Penelope": "LINE/Penelope",
    "PLE": "LINE/Penelope",
    "TIR": "DNA",
    "Helitron": "RC/Helitron",
    "Maverick": "DNA/Maverick",
    "Crypton": "DNA/Crypton",
}

# TEsorter superfamily (or LINE clade) -> RepeatMasker Class/Family, matched
# case-insensitively. Anything unlisted becomes "<RM order>/<TEsorter name>",
# or just the RM order when TEsorter gives no superfamily.
SUPERFAMILY_TO_RM = {
    "gypsy": "LTR/Gypsy",
    "copia": "LTR/Copia",
    "bel-pao": "LTR/Pao",
    "retrovirus": "LTR/ERV",
    "erv": "LTR/ERV",
    "hat": "DNA/hAT",
    "tc1_mariner": "DNA/TcMar",
    "tc1-mariner": "DNA/TcMar",
    "mudr_mutator": "DNA/MULE-MuDR",
    "mutator": "DNA/MULE-MuDR",
    "pif_harbinger": "DNA/PIF-Harbinger",
    "harbinger": "DNA/PIF-Harbinger",
    "piggybac": "DNA/PiggyBac",
    "enspm_cacta": "DNA/CMC-EnSpm",
    "cacta": "DNA/CMC-EnSpm",
    "merlin": "DNA/Merlin",
    "p": "DNA/P",
    "sola": "DNA/Sola",
    "cr1": "LINE/CR1",
    "l1": "LINE/L1",
    "l2": "LINE/L2",
    "rte": "LINE/RTE",
    "i": "LINE/I",
    "jockey": "LINE/I-Jockey",
    "r2": "LINE/R2",
    "rex": "LINE/Rex-Babar",
}

NO_SUPERFAMILY = {"", "unknown", "mixture", "none", "na"}


def bare_name(header):
    return header.split()[0].split("#", 1)[0]


def class_of(header):
    tok = header.split()[0]
    return tok.split("#", 1)[1] if "#" in tok else ""


# ---------------------------------------------------------------- extract

def cmd_extract(args):
    n = n_unknown = n_skipped = 0
    prefixes = tuple(args.keep_prefix)
    classes = open(args.classes_out, "w") if args.classes_out else None
    if classes:
        classes.write("family\trm_class\n")
    with open(args.out, "w") as out:
        for header, seq in iter_fasta(args.library):
            cls = class_of(header)
            if prefixes and not bare_name(header).startswith(prefixes):
                n_skipped += 1
                continue
            if args.unknown_only and cls != "Unknown":
                continue
            write_fasta(out, bare_name(header), seq)
            if classes:
                classes.write(f"{bare_name(header)}\t{cls}\n")
            n += 1
            n_unknown += cls == "Unknown"
    if classes:
        classes.close()
    print(f"[reclassify_unknown] extracted {n} consensi ({n_unknown} Unknown)"
          + (f"; skipped {n_skipped} without a de novo prefix {list(prefixes)} (e.g. Dfam entries)"
             if prefixes else ""), file=sys.stderr)


def read_classes(path):
    """{family: RepeatModeler class} from extract --classes-out."""
    out = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                out[f[0]] = f[1]
    return out


# ---------------------------------------------------------------- readers

def read_tesorter_cls(path):
    """{family: (order, superfamily, clade, n_domains, complete)} from a
    TEsorter *.cls.tsv (#TE Order Superfamily Clade Complete Strand Domains)."""
    calls = {}
    if not path or not os.path.exists(path):
        return calls
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            fam = bare_name(f[0])
            order, superfam, clade = f[1], f[2], f[3]
            complete = f[4] if len(f) > 4 else ""
            domains = [d for d in (f[6].split() if len(f) > 6 else []) if d]
            calls[fam] = (order, superfam, clade, len(domains), complete)
    return calls


def read_diamond(path, max_evalue):
    """{family: (best subject, title, coverage, evidence string)} for every
    query with a host hit (coverage filtering is the caller's). outfmt 6 columns:
    qseqid sseqid pident length qstart qend qlen evalue bitscore stitle.
    Coverage = fraction of the consensus covered by all HSPs to the best
    (highest total bitscore) subject."""
    by_q = {}
    if not path or not os.path.exists(path):
        return {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            q, s = f[0], f[1]
            qs, qe, qlen = int(f[4]), int(f[5]), int(f[6])
            ev, bits = float(f[7]), float(f[8])
            title = f[9] if len(f) > 9 else ""
            if ev > max_evalue:
                continue
            d = by_q.setdefault(bare_name(q), {})
            rec = d.setdefault(s, {"bits": 0.0, "ivs": [], "qlen": qlen, "title": title, "pident": float(f[2])})
            rec["bits"] += bits
            rec["ivs"].append(tuple(sorted((qs, qe))))
    hits = {}
    for q, subs in by_q.items():
        s, rec = max(subs.items(), key=lambda kv: kv[1]["bits"])
        covered = 0
        last_end = 0
        for b, e in sorted(rec["ivs"]):
            b = max(b, last_end + 1)
            if e >= b:
                covered += e - b + 1
            last_end = max(last_end, e)
        cov = covered / rec["qlen"] if rec["qlen"] else 0.0
        hits[q] = (s, rec["title"], cov, f"{s} {cov:.2f}cov {rec['pident']:.0f}%id")
    return hits


def rfam_class(model):
    m = model.lower()
    if "trna" in m:
        return "tRNA"
    if "rrna" in m or m in ("ssu_rrna_eukarya", "lsu_rrna_eukarya"):
        return "rRNA"
    if "srp" in m or "7sk" in m:
        return "srpRNA"
    if re.match(r"^u\d", m) or "snor" in m or "snrna" in m or m.startswith("sno"):
        return "snRNA"
    return "RNA"


def read_rfam(path, cons_len):
    """{family: (rna class, evidence string, coverage)} for every query with an
    included hit (coverage filtering is the caller's), from cmscan --tblout (default
    format): target(model) acc query acc mdl mdl_from mdl_to seq_from seq_to
    strand trunc pass gc bias score E inc description."""
    best = {}
    if not path or not os.path.exists(path):
        return {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            if len(f) < 17:
                continue
            model, q = f[0], bare_name(f[2])
            sf, st = sorted((int(f[7]), int(f[8])))
            score = float(f[14])
            if f[16] != "!":
                continue
            span = st - sf + 1
            if q not in best or score > best[q][0]:
                best[q] = (score, model, span)
    out = {}
    for q, (score, model, span) in best.items():
        L = cons_len.get(q, 0)
        cov = span / L if L else 0.0
        out[q] = (rfam_class(model), f"{model} {cov:.2f}cov {score:.0f}bits", cov)
    return out


def tesorter_to_rm(call, min_domains):
    """RepeatMasker Class/Family for a TEsorter call, or None."""
    if call is None:
        return None
    order, superfam, clade, n_dom, _complete = call
    if n_dom < min_domains:
        return None
    rm_order = ORDER_TO_RM.get(order)
    if rm_order is None:
        return None
    for name in (superfam, clade):
        key = name.strip().lower()
        if key in NO_SUPERFAMILY:
            continue
        if key in SUPERFAMILY_TO_RM:
            mapped = SUPERFAMILY_TO_RM[key]
            if mapped.split("/", 1)[0] == rm_order.split("/", 1)[0]:
                return mapped
        if "/" not in rm_order:
            return f"{rm_order}/{name.strip()}"
    return rm_order


# ---------------------------------------------------------------- merge

# Mavericks/Polintons carry a retroviral-like integrase. GyDB has no
# Maverick models, so it calls that domain LTR (Gypsy / Retroviridae);
# REXdb's Maverick call (which rests on the other Maverick domains too) wins
# instead of a conflict. Ginger (also a Gypsy-like integrase) stays a
# conflict: a REXdb Ginger call alone is weaker evidence.
INTEGRASE_LIKE_DNA = ("DNA/Maverick", "DNA/Polinton")


def domain_evidence(fam, rex, gydb, min_domains):
    """(domain Class/Family or None, conflict text or "") from the two
    TEsorter databases: the more specific call, unless they disagree on order
    (a REXdb Maverick/Polinton call against a GyDB LTR call is not a
    disagreement: see INTEGRASE_LIKE_DNA)."""
    r_call = tesorter_to_rm(rex.get(fam), min_domains)
    g_call = tesorter_to_rm(gydb.get(fam), min_domains)
    if r_call and g_call and r_call.split("/")[0] != g_call.split("/")[0]:
        if r_call.startswith(INTEGRASE_LIKE_DNA) and g_call.split("/")[0] == "LTR":
            return r_call, ""
        return None, f"REXdb {r_call} vs GyDB {g_call}"
    cands = [c for c in (r_call, g_call) if c]
    return (max(cands, key=lambda c: c.count("/")) if cands else None), ""


def multicopy_families(paths, max_copies):
    """Families with more than max_copies .out hits in any sample."""
    if not paths:
        return set()
    return {fam for per_sp in read_family_tandem(paths).values()
            for fam, (_bp, _t, n_hits) in per_sp.items() if n_hits > max_copies}


def cmd_merge(args):
    cons_len = {}
    for header, seq in iter_fasta(args.consensi):
        cons_len[bare_name(header)] = len(seq)
    # With --classes, the consensi are the whole library; only Unknown
    # families get rows.
    families = sorted(cons_len)
    if args.classes:
        rm_class = read_classes(args.classes)
        families = [f for f in families if rm_class.get(f) == "Unknown"]

    rex = read_tesorter_cls(args.tesorter_rexdb)
    gydb = read_tesorter_cls(args.tesorter_gydb)
    host_all = read_diamond(args.diamond, args.max_evalue)
    rfam_all = read_rfam(args.rfam, cons_len)
    # Only qualifying hits can relabel or conflict; every hit is reported.
    host = {q: h for q, h in host_all.items() if h[2] >= args.min_host_cov}
    rfam = {q: r for q, r in rfam_all.items() if r[2] >= args.min_rfam_cov}
    multicopy = multicopy_families(args.family_tandem, args.host_max_copies)

    counts = {}
    with open(args.out, "w") as out:
        out.write("family\told_class\tnew_class\tevidence\ttesorter_rexdb\ttesorter_gydb\t"
                  "diamond_best\trfam_best\tconflict\tnote\n")
        for fam in families:
            domain_class, db_conflict = domain_evidence(fam, rex, gydb, args.min_domains)
            conflict = [db_conflict] if db_conflict else []
            rna = rfam.get(fam)
            hst = host.get(fam)
            note = ""
            if hst and fam in multicopy:
                hst, note = None, "host_annotated_te_orf"
            if domain_class and hst:
                conflict.append(f"domain {domain_class} and host protein {hst[0]}")
            if domain_class and rna:
                conflict.append(f"domain {domain_class} and Rfam {rna[0]}")

            new_class, evidence = "Unknown", ""
            if not conflict:
                if domain_class:
                    new_class, evidence = domain_class, "domain"
                elif rna:
                    new_class, evidence = rna[0], "rfam"
                elif hst:
                    new_class, evidence = "Other/host_gene", "host_protein"
            key = evidence or ("conflict" if conflict else note or "none")
            counts[key] = counts.get(key, 0) + 1

            def fmt(call):
                return "|".join(str(x) for x in call) if call else ""
            out.write(f"{fam}\tUnknown\t{new_class}\t{evidence}\t{fmt(rex.get(fam))}\t"
                      f"{fmt(gydb.get(fam))}\t{host_all[fam][3] if fam in host_all else ''}\t"
                      f"{rfam_all[fam][1] if fam in rfam_all else ''}\t"
                      f"{'; '.join(conflict)}\t{note}\n")
    print(f"[reclassify_unknown] {len(families)} Unknown families: "
          + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())), file=sys.stderr)


RNA_CLASSES = {"rRNA", "tRNA", "snRNA", "srpRNA", "scRNA", "RNA"}
# RepeatMasker orders treated as the same order as TEsorter's call.
ORDER_ALIASES = {"PLE": "LINE"}


def verify_status(rm_cls, dom, rna, host):
    """Status of one classified family. dom: domain Class/Family or None;
    rna: qualifying Rfam (class, ...) or None; host: qualifying host hit or None."""
    rm_order = rm_cls.split("/", 1)[0]
    rm_order = ORDER_ALIASES.get(rm_order, rm_order)
    if rm_order in RNA_CLASSES:
        if rna:
            return "agree_rfam" if rna[0] == rm_order else "rfam_other_rna"
        return "domain_on_rna_label" if dom else "no_evidence"
    if dom:
        dom_order = dom.split("/", 1)[0]
        if rm_order == "Retroposon" and dom_order == "LINE":
            return "retroposon_vs_line"
        if rm_order != dom_order:
            return "disagree_order"
        a, b = rm_cls.lower(), dom.lower()
        if "/" in rm_cls and "/" in dom and (a.startswith(b) or b.startswith(a)):
            return "agree_superfamily"
        return "agree_order"
    if rna:
        return "rfam"
    if host:
        return "host_protein"
    return "no_evidence"


# host_annotated_te_orf (a multi-copy family whose only evidence is a host
# protein) is support for a TE, not a disagreement.
DISAGREEMENT = {"disagree_order", "retroposon_vs_line", "rfam", "host_protein",
                "rfam_other_rna", "domain_on_rna_label", "domain_conflict"}


def read_family_tandem(paths):
    """{sample: {family: (owned_bp, tandem_family, n_hits)}} from family_tandem.tsv."""
    out = {}
    for path in paths:
        with open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            idx = {h: i for i, h in enumerate(header)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                tf = f[idx["tandem_family"]] == "True" if "tandem_family" in idx else False
                out.setdefault(f[idx["sample"]], {})[f[idx["family"]]] = (
                    int(f[idx["owned_bp"]]), tf, int(f[idx["n_hits"]]))
    return out


def cmd_verify(args):
    from summarize_rm import collapse_class

    cons_len = {}
    for header, seq in iter_fasta(args.consensi):
        cons_len[bare_name(header)] = len(seq)
    rm_class = read_classes(args.classes)
    rex = read_tesorter_cls(args.tesorter_rexdb)
    gydb = read_tesorter_cls(args.tesorter_gydb)
    host_all = read_diamond(args.diamond, args.max_evalue)
    rfam_all = read_rfam(args.rfam, cons_len)
    host = {q: h for q, h in host_all.items() if h[2] >= args.min_host_cov}
    rfam = {q: r for q, r in rfam_all.items() if r[2] >= args.min_rfam_cov}
    bp = read_family_tandem(args.family_tandem)
    sample = sorted(bp)
    multicopy = {fam for per_sp in bp.values() for fam, v in per_sp.items() if v[2] > args.host_max_copies}

    rows = []
    for fam, cls in sorted(rm_class.items()):
        if cls in ("Unknown", "") or fam not in cons_len:
            continue
        dom, db_conflict = domain_evidence(fam, rex, gydb, args.min_domains)
        hst = host.get(fam)
        if hst and fam in multicopy:
            hst = None
        status = "domain_conflict" if db_conflict else verify_status(cls, dom, rfam.get(fam), hst)
        if status == "no_evidence" and fam in host and fam in multicopy:
            status = "host_annotated_te_orf"
        owned = {sp: bp[sp].get(fam, (0, False, 0))[0] for sp in sample}
        tandem = any(bp[sp].get(fam, (0, False, 0))[1] for sp in sample)
        rows.append((fam, cls, dom or "", status, owned, tandem, db_conflict,
                     host_all[fam][3] if fam in host_all else "",
                     rfam_all[fam][1] if fam in rfam_all else ""))

    # per sample x collapsed class x status, bp-weighted
    agg = {}
    class_tot = {}
    for fam, cls, _dom, status, owned, _t, *_ in rows:
        top = collapse_class(cls)
        for sp in sample:
            a = agg.setdefault((sp, top, status), [0, 0])
            a[0] += 1 if owned[sp] else 0
            a[1] += owned[sp]
            class_tot[(sp, top)] = class_tot.get((sp, top), 0) + owned[sp]
    with open(args.out, "w") as out:
        out.write("sample\trm_class\tstatus\tn_families_present\towned_bp\tpct_of_class_bp\n")
        for (sp, top, status), (n, b) in sorted(agg.items()):
            if n == 0:
                continue
            tot = class_tot[(sp, top)]
            out.write(f"{sp}\t{top}\t{status}\t{n}\t{b}\t{100.0 * b / tot if tot else 0.0:.2f}\n")

    te_tops = {"DNA", "LINE", "SINE", "LTR", "RC", "Retroposon", "Other"}
    flagged = [r for r in rows if r[3] in DISAGREEMENT
               or (r[5] and collapse_class(r[1]) in te_tops)]
    flagged.sort(key=lambda r: -sum(r[4].values()))
    with open(args.disagreements_out, "w") as out:
        out.write("family\trm_class\tdomain_call\tstatus\ttandem_family\t"
                  + "".join(f"owned_bp_{sp}\t" for sp in sample)
                  + "tesorter_conflict\tdiamond_best\trfam_best\n")
        for fam, cls, dom, status, owned, tandem, conf, hst, rna in flagged:
            out.write(f"{fam}\t{cls}\t{dom}\t{status}\t{tandem}\t"
                      + "".join(f"{owned[sp]}\t" for sp in sample)
                      + f"{conf}\t{hst}\t{rna}\n")

    n_status = {}
    for r in rows:
        n_status[r[3]] = n_status.get(r[3], 0) + 1
    print(f"[reclassify_unknown] verify: {len(rows)} classified families: "
          + ", ".join(f"{k}={v}" for k, v in sorted(n_status.items()))
          + f"; {len(flagged)} flagged", file=sys.stderr)


def reclassified(path):
    """{family: new_class} for families the table relabels (reader for
    summarize_rm.py)."""
    out = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_fam, i_old, i_new = header.index("family"), header.index("old_class"), header.index("new_class")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[i_new] and f[i_new] != f[i_old]:
                out[f[i_fam]] = f[i_new]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("extract")
    e.add_argument("--library", required=True)
    e.add_argument("--out", required=True)
    e.add_argument("--classes-out", default="", help="family -> RepeatModeler class TSV")
    e.add_argument("--unknown-only", action="store_true")
    e.add_argument("--keep-prefix", action="append", default=[],
                   help="only families whose name starts with this (repeatable), e.g. Esto_")
    m = sub.add_parser("merge")
    m.add_argument("--consensi", required=True, help="extract output")
    m.add_argument("--classes", default="", help="extract --classes-out; restricts rows to Unknown families")
    m.add_argument("--tesorter-rexdb", default="")
    m.add_argument("--tesorter-gydb", default="")
    m.add_argument("--diamond", default="", help="DIAMOND blastx outfmt 6 (see read_diamond)")
    m.add_argument("--rfam", default="", help="cmscan --tblout")
    m.add_argument("--min-domains", type=int, default=1)
    m.add_argument("--min-host-cov", type=float, default=0.3)
    m.add_argument("--max-evalue", type=float, default=1e-10)
    m.add_argument("--min-rfam-cov", type=float, default=0.5)
    m.add_argument("--family-tandem", nargs="*", default=[],
                   help="per-sample family_tandem.tsv (n_hits) for --host-max-copies")
    m.add_argument("--host-max-copies", type=int, default=50,
                   help="host_protein only for families with <= this many .out hits in every sample")
    m.add_argument("--out", required=True)
    v = sub.add_parser("verify")
    v.add_argument("--consensi", required=True)
    v.add_argument("--classes", required=True)
    v.add_argument("--tesorter-rexdb", default="")
    v.add_argument("--tesorter-gydb", default="")
    v.add_argument("--diamond", default="")
    v.add_argument("--rfam", default="")
    v.add_argument("--family-tandem", nargs="+", required=True, help="per-sample family_tandem.tsv (owned_bp)")
    v.add_argument("--min-domains", type=int, default=1)
    v.add_argument("--min-host-cov", type=float, default=0.3)
    v.add_argument("--max-evalue", type=float, default=1e-10)
    v.add_argument("--min-rfam-cov", type=float, default=0.5)
    v.add_argument("--host-max-copies", type=int, default=50)
    v.add_argument("--out", required=True, help="class_verification.tsv")
    v.add_argument("--disagreements-out", required=True, help="class_disagreements.tsv")
    args = ap.parse_args()
    {"extract": cmd_extract, "merge": cmd_merge, "verify": cmd_verify}[args.cmd](args)


if __name__ == "__main__":
    main()
