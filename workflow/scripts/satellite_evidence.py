#!/usr/bin/env python3
"""Satellite / rDNA / mito cross-check of the shared library against
TideCluster arrays (trc_crosscheck.py), ribotin rDNA models and MitoHiFi
mitogenomes. Report-only: proposals are written, nothing is relabelled.

Subcommands
  select   write the blastn inputs:
             --db-fasta   candidate library consensi (lib::<family>), every
                          TRC consensus (trc::<sample>:TRC_n), rDNA models
                          (rdna::<name>) and mitogenomes (mito::<name>);
                          searched all-vs-all
             --refs-fasta rDNA + mito only; the whole library is searched
                          against it (rDNA / NUMT pieces under TE labels)
           Candidates: tandem_family in any sample, >= --min-trc-bp inside
           TRCs in any sample, or listed in --focus.
  report   merge everything into four tables (see README "Satellite, rDNA
           and mito cross-check"):
             satellite_family_evidence.tsv  family x sample
             satellite_family_calls.tsv     one row per family: proposed
                                            class, confidence, reason
             satellite_family_pairs.tsv     same satellite or not?
             satellite_proposals.tsv        curated_families format
Stdlib only."""

import argparse
import sys
from collections import defaultdict

from family_groups import merged_len
from fasta_utils import iter_fasta, seq_id, write_fasta

SEP = "::"


def bare(name):
    return seq_id(name).split("#", 1)[0]


def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if line.strip():
                yield dict(zip(header, line.rstrip("\n").split("\t")))


def num(x, cast=float):
    try:
        return cast(x)
    except (TypeError, ValueError):
        return None


def label_fastas(specs):
    """'label=path' -> [(label, path)]; label is rdna or mito."""
    out = []
    for spec in specs or []:
        label, path = spec.split("=", 1)
        out.append((label, path))
    return out


# ---------------------------------------------------------------- select
def cmd_select(args):
    focus = set(args.focus or [])
    cand = set(focus)
    for r in read_tsv(args.family_tandem):
        if r.get("tandem_family") == "True":
            cand.add(r["family"])
    for path in args.crosscheck or []:
        for r in read_tsv(path):
            if num(r["in_trc_bp"], int) and int(r["in_trc_bp"]) >= args.min_trc_bp:
                cand.add(r["family"])
    n = {"lib": 0, "trc": 0, "rdna": 0, "mito": 0}
    seen_refs = set()
    with open(args.db_fasta, "w") as db, open(args.refs_fasta, "w") as refs:
        for header, seq in iter_fasta(args.library):
            if bare(header) in cand:
                write_fasta(db, f"lib{SEP}{bare(header)}", seq)
                n["lib"] += 1
        for path in args.trc_consensus or []:
            for header, seq in iter_fasta(path):
                write_fasta(db, f"trc{SEP}{seq_id(header)}", seq)
                n["trc"] += 1
        for label, path in label_fastas(args.ref):
            for header, seq in iter_fasta(path):
                name = f"{label}{SEP}{seq_id(header)}"
                if name in seen_refs:  # same model listed for several samples
                    continue
                seen_refs.add(name)
                write_fasta(db, name, seq)
                write_fasta(refs, name, seq)
                n[label] += 1
    missing = focus - {bare(h) for h, _ in iter_fasta(args.library)}
    if missing:
        print(f"[satellite_evidence] WARNING: focus families not in the library: {sorted(missing)}",
              file=sys.stderr)
    print(f"[satellite_evidence] blast db: {n}; {len(cand)} candidate families", file=sys.stderr)


# ---------------------------------------------------------------- report
def read_blast(path):
    """Yield (q, s, pident, qstart, qend, qlen) with the label prefixes kept."""
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            q, s = f[0], f[1]
            yield q, s, float(f[2]), int(f[4]), int(f[5]), int(f[8]), int(f[6]), int(f[7]), int(f[9])


def coverage(hits, min_id):
    """{(query, subject): (query coverage, bp-weighted identity)} over HSPs
    >= min_id; coverage = merged query bp / query length."""
    ivs, idw, qlen = defaultdict(list), defaultdict(float), {}
    for q, s, pid, qs, qe, ql, ss, se, sl in hits:
        if q == s or pid < min_id:
            continue
        ivs[(q, s)].append((min(qs, qe), max(qs, qe)))
        idw[(q, s)] += pid * abs(qe - qs)
        qlen[(q, s)] = ql
    out = {}
    for k, iv in ivs.items():
        span = sum(e - s for s, e in iv) or 1
        out[k] = (merged_len(iv) / qlen[k], idw[k] / span)
    return out


def strip(name):
    return name.split(SEP, 1)[1] if SEP in name else name


def monomer_match(a, b, tol):
    if not a or not b:
        return "NA"
    lo, hi = sorted((a, b))
    if abs(hi - lo) <= tol * hi:
        return "same"
    k = round(hi / lo)
    if k >= 2 and abs(hi - k * lo) <= tol * hi:
        return f"multiple_x{k}"
    return "different"


def cmd_report(args):
    t = args
    ft = defaultdict(dict)  # family -> sample -> family_tandem row
    for r in read_tsv(t.family_tandem):
        if r.get("arm", "shared") == "shared":
            ft[r["family"]][r["species"]] = r
    xc = defaultdict(dict)  # family -> sample -> crosscheck row
    tc_samples = []
    for path in t.crosscheck:
        for r in read_tsv(path):
            xc[r["family"]][r["species"]] = r
            if r["species"] not in tc_samples:
                tc_samples.append(r["species"])
    trc_info = {}
    for path in t.trc_info or []:
        for r in read_tsv(path):
            trc_info[(r["species"], r["trc"])] = r
    max_period = {}
    for path in t.tc_params or []:
        for r in read_tsv(path):
            max_period[r["species"]] = int(r["tidehunter_max_period"])
    infer = dict(x.split("=", 1) for x in (t.infer or []))  # target=donor
    focus = list(t.focus or [])
    expected_independent = {frozenset(x.split(",")) for x in (t.expected_independent or [])}

    sat_hits = list(read_blast(t.sat_blast))
    ref_hits = list(read_blast(t.ref_blast)) if t.ref_blast else []
    cov_lib = coverage(sat_hits, 0.0)  # all-vs-all, identity reported
    cov_rdna = coverage([h for h in ref_hits if h[1].startswith("rdna" + SEP)], t.rdna_min_id)
    cov_mito = coverage([h for h in ref_hits if h[1].startswith("mito" + SEP)], t.mito_min_id)

    def best(cov, family, prefix):
        hits = [(c, i, s) for (q, s), (c, i) in cov.items()
                if bare(q.split(SEP, 1)[-1]) == family and s.startswith(prefix + SEP)]
        return max(hits, default=None)

    families = sorted(set(xc) | {f for f, d in ft.items() if any(r.get("tandem_family") == "True" for r in d.values())}
                      | set(focus)
                      | {bare(q) for (q, s) in list(cov_rdna) + list(cov_mito)})

    # ---- family x sample evidence
    ev_rows = []
    for fam in families:
        cf = next((r["class_family"] for r in ft.get(fam, {}).values()), None) or \
             next((r["class_family"] for r in xc.get(fam, {}).values()), "NA")
        for s in tc_samples + [x for x in infer if x not in tc_samples]:
            src = "tidecluster" if s in tc_samples else f"inferred:{infer[s]}"
            ftr = ft.get(fam, {}).get(s, {})
            x = xc.get(fam, {}).get(s if s in tc_samples else infer[s], {})
            ti = trc_info.get((s if s in tc_samples else infer[s], x.get("top_trc")), {})
            ev_rows.append({
                "family": fam, "class_family": cf, "sample": s, "evidence_source": src,
                "owned_bp": ftr.get("owned_bp", "0"), "tandem_family": ftr.get("tandem_family", "NA"),
                "family_tandem_monomer_period": ftr.get("monomer_period", "NA"),
                "frac_in_trc": x.get("frac_in_trc", "0.000" if s in tc_samples else "NA"),
                "top_trc": x.get("top_trc", "NA"), "top_trc_bp": x.get("top_trc_bp", "0"),
                "trc_superfamily": ti.get("superfamily", "NA"), "trc_rdna_flag": ti.get("rdna_flag", "NA"),
                "tarean_monomer_len": ti.get("tarean_monomer_len", "NA"),
                "th_monomer_median": x.get("th_monomer_median", "NA"),
                "th_monomer_support": x.get("th_monomer_support", "NA"),
            })
    ev_cols = list(ev_rows[0].keys()) if ev_rows else ["family"]
    with open(t.evidence_out, "w") as out:
        out.write("\t".join(ev_cols) + "\n")
        for r in ev_rows:
            out.write("\t".join(str(r[c]) for c in ev_cols) + "\n")
    ev = defaultdict(dict)
    for r in ev_rows:
        ev[r["family"]][r["sample"]] = r

    # ---- per-family calls
    calls = {}
    with open(t.calls_out, "w") as out:
        out.write("family\tclass_family\tsat_samples_pass\tsat_samples_eligible\tbeyond_tidehunter\tinferred_support\t"
                  "best_trc_hit\trdna_cov\trdna_id\tmito_cov\tmito_id\tproposed_class_family\t"
                  "confidence\treason\n")
        for fam in families:
            rows = ev[fam]
            cf = next(iter(rows.values()))["class_family"] if rows else "NA"
            elig, passed, weak, beyond = [], [], [], []
            for s in tc_samples:
                r = rows.get(s)
                if not r:
                    continue
                owned, frac = num(r["owned_bp"], int) or 0, num(r["frac_in_trc"]) or 0.0
                ok = r["tandem_family"] == "True" and frac >= t.min_trc_cov
                # A tandem family whose repeat unit is longer than TideHunter looked
                # for can't be in a TRC: no TideCluster evidence either way.
                ftr = ft.get(fam, {}).get(s, {})
                unit = num(ftr.get("monomer_period")) or num(ftr.get("cons_len")) or 0
                too_long = (not ok and r["tandem_family"] == "True"
                            and (num(ftr.get("tandem_frac")) or 0) >= t.min_tandem_frac_long
                            and unit > max_period.get(s, 3000))
                if owned >= t.major_min_bp:
                    elig.append(s)
                    if ok:
                        passed.append(s)
                    elif too_long:
                        beyond.append(s)
                elif ok or frac >= 0.3:
                    weak.append(s)
            inferred = []
            for target, donor in infer.items():
                r = rows.get(target)
                if not r:
                    continue
                if donor in passed:
                    has = r["tandem_family"] == "True" and (num(r["owned_bp"], int) or 0) >= t.major_min_bp
                    inferred.append(f"{target}:{'confirms' if has else 'no_arrays'}")
            b_trc = best(cov_lib, fam, "trc")
            b_rd = best(cov_rdna, fam, "rdna")
            b_mt = best(cov_mito, fam, "mito")
            rd_flag = any(r["trc_rdna_flag"] == "True" for r in rows.values())

            prop, conf, why = "NA", "none", []
            if b_rd and b_rd[0] >= t.rdna_min_cov:
                prop, conf = "rRNA", "high"
                why.append(f"{b_rd[0]:.2f} of consensus matches ribotin {strip(b_rd[2])} at {b_rd[1]:.1f}%")
            elif b_mt and b_mt[0] >= t.mito_min_cov:
                prop, conf = "Other/NUMT", "high"
                why.append(f"{b_mt[0]:.2f} of consensus matches mitogenome {strip(b_mt[2])} at {b_mt[1]:.1f}%")
            elif elig and passed and len(passed) == len(elig):
                prop, conf = "Satellite", "high"
                why.append(f"tandem and >= {t.min_trc_cov} of its bp in TRC arrays in all eligible "
                           f"samples ({','.join(passed)})")
            elif elig and len(passed) + len(beyond) == len(elig):
                prop, conf = "Satellite", "medium"
                why.append(f"repeat unit longer than TideHunter's max period in {','.join(beyond)} "
                           f"(no TRC possible); tandem by RepeatMasker arrays (family_tandem) there"
                           + (f"; TRC support in {','.join(passed)}" if passed else ""))
            elif passed or weak:
                prop, conf = "Satellite", "medium" if passed else "low"
                why.append(f"TRC support in {','.join(passed + weak)} of eligible {','.join(elig) or 'none'}")
            if prop in ("Satellite",) and rd_flag:
                prop, conf = "rRNA", "medium"
                why.append("its top TRC is flagged rDNA by TideCluster (no ribotin match)")
            if prop != "NA" and cf not in ("NA", "Unknown", prop) and not cf.startswith(prop):
                why.append(f"currently {cf} ({cf}-derived)")
            calls[fam] = {"class_family": cf, "prop": prop, "conf": conf}
            out.write(f"{fam}\t{cf}\t{len(passed)}\t{len(elig)}\t{','.join(beyond) or '.'}\t{','.join(inferred) or '.'}\t"
                      + (f"{strip(b_trc[2])}:cov={b_trc[0]:.2f}:id={b_trc[1]:.1f}" if b_trc else "NA") + "\t"
                      + (f"{b_rd[0]:.2f}\t{b_rd[1]:.1f}" if b_rd else "NA\tNA") + "\t"
                      + (f"{b_mt[0]:.2f}\t{b_mt[1]:.1f}" if b_mt else "NA\tNA") + "\t"
                      + f"{prop}\t{conf}\t{'; '.join(why) or '.'}\n")

    # ---- pairs
    def trc_bp(fam, s):
        x = xc.get(fam, {}).get(s)
        if not x or x["trc_bp"] == ".":
            return {}
        return {k: int(v) for k, v in (kv.rsplit(":", 1) for kv in x["trc_bp"].split(";"))}

    def monomer(fam):
        best_s = max(tc_samples, key=lambda s: num(xc.get(fam, {}).get(s, {}).get("tidehunter_bp"), int) or 0,
                     default=None)
        m = num(xc.get(fam, {}).get(best_s, {}).get("th_monomer_median")) if best_s else None
        if m:
            return m, f"tidehunter:{best_s}"
        p = next((num(r.get("monomer_period")) for r in ft.get(fam, {}).values()
                  if num(r.get("monomer_period"))), None)
        return (p, "consensus_period") if p else (None, "NA")

    def superfam(fam):
        out = set()
        for s in tc_samples:
            r = ev[fam].get(s)
            if r and r["trc_superfamily"] != "NA":
                out.add((s, r["trc_superfamily"]))
        return out

    tandemish = {f for f in families
                 if (calls[f]["prop"] in ("Satellite", "rRNA") and calls[f]["conf"] in ("high", "medium"))
                 or f in focus}
    pairs = set()
    for i, a in enumerate(focus):
        for b in focus[i + 1:]:
            pairs.add(tuple(sorted((a, b))))
    for s in tc_samples:
        by_trc = defaultdict(set)
        for f in tandemish:
            for trc in trc_bp(f, s):
                by_trc[trc].add(f)
        for fs in by_trc.values():
            fs = sorted(fs)
            for i, a in enumerate(fs):
                for b in fs[i + 1:]:
                    pairs.add((a, b))
    for (q, s) in cov_lib:
        if q.startswith("lib" + SEP) and s.startswith("lib" + SEP):
            a, b = strip(q), strip(s)
            if a in tandemish and b in tandemish:
                pairs.add(tuple(sorted((a, b))))

    groups = {}  # union-find over same_satellite pairs, for the proposals

    def find(x):
        while groups.get(x, x) != x:
            x = groups[x]
        return x

    with open(t.pairs_out, "w") as out:
        out.write("family_a\tfamily_b\tcov_a_by_b\tcov_b_by_a\tidentity\tshared_trc_frac\t"
                  "max_shared_trc_frac\tmonomer_a\tmonomer_b\tmonomer_match\tsame_trc_superfamily\t"
                  "verdict\tnote\n")
        for a, b in sorted(pairs):
            ab = cov_lib.get((f"lib{SEP}{a}", f"lib{SEP}{b}"))
            ba = cov_lib.get((f"lib{SEP}{b}", f"lib{SEP}{a}"))
            cov_ab, cov_ba = (ab[0] if ab else 0.0), (ba[0] if ba else 0.0)
            ident = max([x[1] for x in (ab, ba) if x], default=0.0)
            shared, per = [], []
            for s in tc_samples:
                ta, tb = trc_bp(a, s), trc_bp(b, s)
                if not ta or not tb:
                    per.append(f"{s}:NA")
                    continue
                v = sum(min(ta[k], tb[k]) for k in set(ta) & set(tb)) / min(sum(ta.values()), sum(tb.values()))
                shared.append(v)
                per.append(f"{s}:{v:.2f}")
            mx = max(shared, default=None)
            (ma, sa), (mb, sb) = monomer(a), monomer(b)
            mm = monomer_match(ma, mb, t.monomer_tol)
            sfa, sfb = superfam(a), superfam(b)
            same_sf = "NA" if not (sfa and sfb) else str(bool(sfa & sfb))
            cov = max(cov_ab, cov_ba)
            sim = ident >= t.min_pair_id and cov >= t.min_pair_cov
            related = ident >= t.related_min_id and cov >= t.min_pair_cov
            co = mx is not None and mx >= t.min_shared_trc_frac
            mono_ok = mm == "same" or mm.startswith("multiple")
            if co and sim and mono_ok:
                verdict = "same_satellite"
            elif co and related and mono_ok:
                verdict = "same_satellite_diverged"  # same arrays + monomer, consensi < min_pair_id
            elif co and related:
                verdict = "co_located_related" if mm == "different" else "undetermined"
            elif co:
                verdict = "co_located_distinct"
            elif related and mx is not None:
                verdict = "related_not_co_located"
            elif related:
                verdict = "related_no_shared_sample"
            elif mx is not None and mx <= t.max_independent_shared_frac and cov < t.independent_max_cov:
                verdict = "independent"
            else:
                verdict = "undetermined"
            note = []
            if frozenset((a, b)) in expected_independent and verdict == "same_satellite":
                note.append("CONTRADICTS expected_independent: check by hand before grouping")
            elif frozenset((a, b)) in expected_independent:
                note.append("expected_independent")
            if mm.startswith("multiple"):
                note.append("monomers differ by an integer factor (HOR-like)")
            if verdict == "same_satellite" and not note:
                ra, rb = find(a), find(b)
                if ra != rb:
                    groups[rb] = ra
            out.write(f"{a}\t{b}\t{cov_ab:.2f}\t{cov_ba:.2f}\t{ident:.1f}\t{';'.join(per)}\t"
                      f"{'NA' if mx is None else f'{mx:.2f}'}\t"
                      f"{'NA' if ma is None else f'{ma:g}'}({sa})\t{'NA' if mb is None else f'{mb:g}'}({sb})\t"
                      f"{mm}\t{same_sf}\t{verdict}\t{'; '.join(note) or '.'}\n")

    # ---- proposals (curated_families format)
    with open(t.proposals_out, "w") as out:
        out.write("group\tfamily\tclass_family\tpart\tevidence\tnote\tconfidence\n")
        for fam in families:
            c = calls[fam]
            if c["conf"] not in ("high", "medium"):
                continue
            # same_satellite pairs share one group name (never merged, only grouped)
            group = f"{c['prop'].split('/')[-1]}-{find(fam)}"
            note = f"{c['class_family']}-derived" if c["class_family"] not in ("NA", "Unknown", c["prop"]) else "."
            out.write(f"{group}\t{fam}\t{c['prop']}\t"
                      f"{'45S_unit' if c['prop'] == 'rRNA' else '.'}\tsatellite_crosscheck\t{note}\t{c['conf']}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("select")
    s.add_argument("--library", required=True)
    s.add_argument("--family-tandem", required=True, help="combined shared-arm family_tandem.tsv")
    s.add_argument("--crosscheck", nargs="*", default=[], help="trc_crosscheck.py --family-out tables")
    s.add_argument("--trc-consensus", nargs="*", default=[])
    s.add_argument("--ref", nargs="*", default=[], help="rdna=path / mito=path")
    s.add_argument("--focus", nargs="*", default=[])
    s.add_argument("--min-trc-bp", type=int, default=10000)
    s.add_argument("--db-fasta", required=True)
    s.add_argument("--refs-fasta", required=True)

    r = sub.add_parser("report")
    r.add_argument("--family-tandem", required=True)
    r.add_argument("--crosscheck", nargs="+", required=True)
    r.add_argument("--trc-info", nargs="*", default=[])
    r.add_argument("--sat-blast", required=True)
    r.add_argument("--ref-blast", default="")
    r.add_argument("--infer", nargs="*", default=[], help="target=donor sample pairs, e.g. ind8=ind6")
    r.add_argument("--focus", nargs="*", default=[])
    r.add_argument("--expected-independent", nargs="*", default=[], help="famA,famB")
    r.add_argument("--min-trc-cov", type=float, default=0.7)
    r.add_argument("--major-min-bp", type=int, default=100000)
    r.add_argument("--min-pair-id", type=float, default=80.0)
    r.add_argument("--min-pair-cov", type=float, default=0.5)
    r.add_argument("--related-min-id", type=float, default=65.0)
    r.add_argument("--independent-max-cov", type=float, default=0.2)
    r.add_argument("--min-tandem-frac-long", type=float, default=0.9)
    r.add_argument("--tc-params", nargs="*", default=[], help="tidecluster_regions.py --params tables")
    r.add_argument("--min-shared-trc-frac", type=float, default=0.5)
    r.add_argument("--max-independent-shared-frac", type=float, default=0.1)
    r.add_argument("--monomer-tol", type=float, default=0.05)
    r.add_argument("--rdna-min-cov", type=float, default=0.5)
    r.add_argument("--rdna-min-id", type=float, default=90.0)
    r.add_argument("--mito-min-cov", type=float, default=0.5)
    r.add_argument("--mito-min-id", type=float, default=80.0)
    r.add_argument("--evidence-out", required=True)
    r.add_argument("--calls-out", required=True)
    r.add_argument("--pairs-out", required=True)
    r.add_argument("--proposals-out", required=True)

    args = ap.parse_args()
    cmd_select(args) if args.cmd == "select" else cmd_report(args)


if __name__ == "__main__":
    main()
