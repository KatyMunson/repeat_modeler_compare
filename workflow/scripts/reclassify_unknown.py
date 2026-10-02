#!/usr/bin/env python3
"""Reclassify the library's Unknown families from three cheap, independent
screens of their consensus sequences:

  TEsorter   TE protein domains (REXdb metazoa, GyDB) -> order / superfamily
  Rfam       structured RNA (cmscan --cut_ga)         -> rRNA / tRNA / snRNA / ...
  DIAMOND    host proteins (TE proteins removed)      -> Other/host_gene

Subcommands:
  extract   write every #Unknown library record, header = bare family name
            (TEsorter and cmscan mangle '#'), to a FASTA
  merge     combine the screens into one reclassification table

Precedence (merge), first match wins:
  1. TEsorter call at order level or finer, REXdb and GyDB not disagreeing
     on order, and >= --min-domains domains           evidence = domain
  2. Rfam hit covering >= --min-rfam-cov of the consensus  evidence = rfam
  3. DIAMOND host hit covering >= --min-host-cov       evidence = host_protein
  Disagreements (REXdb vs GyDB order; a domain call AND a host hit; a domain
  call AND an Rfam hit) leave the family Unknown with `conflict` recorded.
A domain call beats Rfam and host hits only when they don't also qualify:
both qualifying is a conflict, not a win, because a host gene or RNA
carrying a TE domain (a domesticated TE, a TE-derived exon) is ambiguous.

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
    n = 0
    with open(args.out, "w") as out:
        for header, seq in iter_fasta(args.library):
            if class_of(header) != "Unknown":
                continue
            write_fasta(out, bare_name(header), seq)
            n += 1
    print(f"[reclassify_unknown] extracted {n} Unknown consensi", file=sys.stderr)


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

def cmd_merge(args):
    cons_len = {}
    for header, seq in iter_fasta(args.consensi):
        cons_len[bare_name(header)] = len(seq)

    rex = read_tesorter_cls(args.tesorter_rexdb)
    gydb = read_tesorter_cls(args.tesorter_gydb)
    host_all = read_diamond(args.diamond, args.max_evalue)
    rfam_all = read_rfam(args.rfam, cons_len)
    # Only qualifying hits can relabel or conflict; every hit is reported.
    host = {q: h for q, h in host_all.items() if h[2] >= args.min_host_cov}
    rfam = {q: r for q, r in rfam_all.items() if r[2] >= args.min_rfam_cov}

    counts = {}
    with open(args.out, "w") as out:
        out.write("family\told_class\tnew_class\tevidence\ttesorter_rexdb\ttesorter_gydb\t"
                  "diamond_best\trfam_best\tconflict\n")
        for fam in sorted(cons_len):
            r_call = tesorter_to_rm(rex.get(fam), args.min_domains)
            g_call = tesorter_to_rm(gydb.get(fam), args.min_domains)
            conflict = []
            domain_class = None
            if r_call and g_call and r_call.split("/")[0] != g_call.split("/")[0]:
                conflict.append(f"REXdb {r_call} vs GyDB {g_call}")
            else:
                # prefer the more specific call (REXdb first when equal)
                cands = [c for c in (r_call, g_call) if c]
                if cands:
                    domain_class = max(cands, key=lambda c: c.count("/"))
            rna = rfam.get(fam)
            hst = host.get(fam)
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
            counts[evidence or ("conflict" if conflict else "none")] = \
                counts.get(evidence or ("conflict" if conflict else "none"), 0) + 1

            def fmt(call):
                return "|".join(str(x) for x in call) if call else ""
            out.write(f"{fam}\tUnknown\t{new_class}\t{evidence}\t{fmt(rex.get(fam))}\t"
                      f"{fmt(gydb.get(fam))}\t{host_all[fam][3] if fam in host_all else ''}\t"
                      f"{rfam_all[fam][1] if fam in rfam_all else ''}\t"
                      f"{'; '.join(conflict)}\n")
    print(f"[reclassify_unknown] {len(cons_len)} Unknown families: "
          + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())), file=sys.stderr)


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
    m = sub.add_parser("merge")
    m.add_argument("--consensi", required=True, help="extract output")
    m.add_argument("--tesorter-rexdb", default="")
    m.add_argument("--tesorter-gydb", default="")
    m.add_argument("--diamond", default="", help="DIAMOND blastx outfmt 6 (see read_diamond)")
    m.add_argument("--rfam", default="", help="cmscan --tblout")
    m.add_argument("--min-domains", type=int, default=1)
    m.add_argument("--min-host-cov", type=float, default=0.3)
    m.add_argument("--max-evalue", type=float, default=1e-10)
    m.add_argument("--min-rfam-cov", type=float, default=0.5)
    m.add_argument("--out", required=True)
    args = ap.parse_args()
    {"extract": cmd_extract, "merge": cmd_merge}[args.cmd](args)


if __name__ == "__main__":
    main()
