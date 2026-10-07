#!/usr/bin/env python3
"""Satellite / rDNA / mito cross-check of the shared library against
TideCluster arrays (trc_crosscheck.py), ribotin rDNA models, MitoHiFi
mitogenomes and, optionally, a harmonized satellite motif library (the
satellite pipeline's harmonized_repeatmasker_lib.fasta). Report-only:
proposals are written, nothing is relabelled.

Two kinds of satellite evidence, used together:
  arrays   TideCluster: is the family's sequence in tandem arrays *in this
           assembly* (frac_in_trc)? Needs only a TideCluster run.
  motif    harmonized library: is the consensus made of a known satellite
           unit? Names the satellite (group = motif ID) and catches
           families whose arrays TideCluster didn't cluster, whose copies
           aren't flagged tandem_family, or that are below major_min_bp.

Subcommands
  select   write the blastn inputs:
             --db-fasta   candidate library consensi (lib::<family>), every
                          TRC consensus (trc::<sample>:TRC_n), rDNA models
                          (rdna::<name>), mitogenomes (mito::<name>) and,
                          with --motifs, every harmonized motif tiled to
                          >= --motif-tile bp (motif::<name>, rotation-safe
                          like the TRC dimers); searched all-vs-all
             --refs-fasta rDNA + mito only; the whole library is searched
                          against it (rDNA / NUMT pieces under TE labels)
           Candidates: tandem_family in any sample, >= --min-trc-bp inside
           TRCs in any sample, or listed in --focus.
  apply    satellite_proposals.tsv -> the rows summarize applies
           (satellite_crosscheck.apply): only --confidence levels, minus
           families the user's curated table already labels (theirs wins)
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
TE_CLASSES = ("DNA", "LINE", "SINE", "LTR", "RC", "Retroposon", "PLE")


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
    unit = {}  # family -> repeat unit (KITE founder > consensus self-period)
    for r in read_tsv(args.family_tandem):
        if r.get("tandem_family") == "True":
            cand.add(r["family"])
        if num(r.get("monomer_period")):
            unit.setdefault(r["family"], num(r["monomer_period"]))
    for path in args.crosscheck or []:
        for r in read_tsv(path):
            if num(r["in_trc_bp"], int) and int(r["in_trc_bp"]) >= args.min_trc_bp:
                cand.add(r["family"])
            if num(r.get("kite_founder_median")):
                unit[r["family"]] = num(r["kite_founder_median"])
    n = {"lib": 0, "trc": 0, "rdna": 0, "mito": 0}
    seen_refs = set()
    with open(args.db_fasta, "w") as db, open(args.refs_fasta, "w") as refs:
        for header, seq in iter_fasta(args.library):
            fam = bare(header)
            if fam in cand:
                # About one monomer long: write it as a dimer so any rotation of
                # the other sequence aligns in one piece (folded back in report).
                if unit.get(fam) and len(seq) < 1.5 * unit[fam]:
                    write_fasta(db, f"libx2{SEP}{fam}", seq + seq)
                    n["dimer"] = n.get("dimer", 0) + 1
                    print(f"[satellite_evidence] {fam}: {len(seq)} bp ~ one {unit[fam]:g} bp unit; "
                          f"searched as a dimer", file=sys.stderr)
                else:
                    write_fasta(db, f"lib{SEP}{fam}", seq)
                n["lib"] += 1
        for path in args.trc_consensus or []:
            for header, seq in iter_fasta(path):
                write_fasta(db, f"trc{SEP}{seq_id(header)}", seq)
                n["trc"] += 1
        if args.motifs:
            for header, seq in iter_fasta(args.motifs):
                if not seq:
                    continue
                reps = max(2, -(-args.motif_tile // len(seq)))
                write_fasta(db, f"motif{SEP}{bare(header)}", seq * reps)
                n["motif"] = n.get("motif", 0) + 1
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
    """Yield (q, s, pident, q_intervals, qlen) with the label prefixes kept.
    Library consensi searched as dimers (libx2::) are renamed lib:: and their
    query coordinates folded back onto the original length."""
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            q, s = f[0], f[1].replace(f"libx2{SEP}", f"lib{SEP}", 1)
            qs, qe = sorted((int(f[4]), int(f[5])))
            ql = int(f[8])
            ivs = [(qs, qe)]
            if q.startswith(f"libx2{SEP}"):
                q = q.replace(f"libx2{SEP}", f"lib{SEP}", 1)
                ql //= 2
                ivs = []
                for a, b in [(qs, min(qe, ql)), (max(qs, ql + 1) - ql, qe - ql)]:
                    if a <= b:
                        ivs.append((max(a, 1), min(b, ql)))
                if qe - qs + 1 >= ql:
                    ivs = [(1, ql)]
            yield q, s, float(f[2]), ivs, ql


def coverage(hits, min_id):
    """{(query, subject): (query coverage, bp-weighted identity)} over HSPs
    >= min_id; coverage = merged query bp / query length."""
    ivs, idw, qlen = defaultdict(list), defaultdict(float), {}
    for q, s, pid, q_ivs, ql in hits:
        if q == s or pid < min_id:
            continue
        for a, b in q_ivs:
            ivs[(q, s)].append((a, b))
            idw[(q, s)] += pid * (b - a)
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
    # HOR-like: hi is k copies of lo, within tol of the *shorter* unit (a
    # tolerance on hi would make any long period a "multiple" at large k)
    k = round(hi / lo)
    if k >= 2 and abs(hi - k * lo) <= tol * lo:
        return f"multiple_x{k}"
    return "different"


def cmd_report(args):
    t = args
    ft = defaultdict(dict)  # family -> sample -> family_tandem row
    for r in read_tsv(t.family_tandem):
        if r.get("arm", "shared") == "shared":
            ft[r["family"]][r["sample"]] = r
    xc = defaultdict(dict)  # family -> sample -> crosscheck row
    tc_samples = []
    for path in t.crosscheck:
        for r in read_tsv(path):
            xc[r["family"]][r["sample"]] = r
            if r["sample"] not in tc_samples:
                tc_samples.append(r["sample"])
    trc_info = {}
    for path in t.trc_info or []:
        for r in read_tsv(path):
            trc_info[(r["sample"], r["trc"])] = r
    max_period = {}
    for path in t.tc_params or []:
        for r in read_tsv(path):
            max_period[r["sample"]] = int(r["tidehunter_max_period"])
    infer = dict(x.split("=", 1) for x in (t.infer or []))  # target=donor
    scale = {k: int(v) for k, v in (x.split("=", 1) for x in (t.copy_scale or []))}
    focus = list(t.focus or [])
    expected_independent = {frozenset(x.split(",")) for x in (t.expected_independent or [])}

    sat_hits = list(read_blast(t.sat_blast))
    ref_hits = list(read_blast(t.ref_blast)) if t.ref_blast else []
    cov_lib = coverage(sat_hits, 0.0)  # all-vs-all, identity reported
    cov_rdna = coverage([h for h in ref_hits if h[1].startswith("rdna" + SEP)], t.rdna_min_id)
    cov_mito = coverage([h for h in ref_hits if h[1].startswith("mito" + SEP)], t.mito_min_id)
    # library family covered by tiled motifs: HSPs >= motif_min_id
    cov_motif = coverage([h for h in sat_hits if h[0].startswith("lib" + SEP) and h[1].startswith("motif" + SEP)],
                         t.motif_min_id)
    # Best motif per family: coverage weighted by identity above the floor,
    # so a near-identical long unit over half the consensus (MLI_SAT728_a,
    # 0.45 at 99.7%) beats a short motif loosely tiled over more of it
    # (MLI_SAT28_ap, 0.58 at 81%); plain cov * identity picked the latter.
    def motif_score(c, i):
        return c * (i - t.motif_min_id + 1)
    motif_best = {}
    for (q, sub), (c, i) in cov_motif.items():
        fam = bare(strip(q))
        if fam not in motif_best or motif_score(c, i) > motif_score(*motif_best[fam][:2]):
            motif_best[fam] = (c, i, strip(sub))

    def motif_tier(fam):
        m = motif_best.get(fam)
        if not m:
            return None, m
        c, i, _ = m
        if c >= t.motif_high_cov and i >= t.motif_high_id:
            return "high", m
        if c >= t.motif_min_cov:
            return "medium", m
        if c >= t.motif_partial_cov and i >= t.motif_high_id:
            return "partial", m
        return None, m

    def best(cov, family, prefix):
        hits = [(c, i, s) for (q, s), (c, i) in cov.items()
                if bare(q.split(SEP, 1)[-1]) == family and s.startswith(prefix + SEP)]
        return max(hits, default=None)

    families = sorted(set(xc) | {f for f, d in ft.items() if any(r.get("tandem_family") == "True" for r in d.values())}
                      | set(focus)
                      | {bare(q) for (q, s) in list(cov_rdna) + list(cov_mito)}
                      | set(motif_best))

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
                "tidehunter_bp": x.get("tidehunter_bp", "0"),
                "th_monomer_median": x.get("th_monomer_median", "NA"),
                "th_monomer_support": x.get("th_monomer_support", "NA"),
                "kite_founder": x.get("kite_founder_median", "NA"),
                "kite_founder_support": x.get("kite_founder_support", "NA"),
                "trc_kite_founder": ti.get("kite_founder_median", "NA"),
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
                  "best_trc_hit\tbest_motif\tmotif_cov\tmotif_id\tmotif_tier\trdna_cov\trdna_id\tmito_cov\tmito_id\t"
                  "proposed_class_family\t"
                  "confidence\treason\n")
        for fam in families:
            rows = ev[fam]
            cf = next(iter(rows.values()))["class_family"] if rows else "NA"
            elig, passed, weak, beyond, partial = [], [], [], [], []
            # Repeat unit for the TideHunter-limit test: KITE founder (any
            # TideCluster sample) > consensus self-period > consensus length.
            kite_unit = next((num(xc.get(fam, {}).get(s, {}).get("kite_founder_median")) for s in tc_samples
                              if num(xc.get(fam, {}).get(s, {}).get("kite_founder_median"))), None)
            # A sample only counts toward "all eligible samples" if the family
            # is tandem there or holds a real share of its bp there: a few
            # scattered copies in another species must not veto a satellite
            # that fills arrays in one (e.g. TRC_1: 542 Mb vs 0.1 Mb).
            top_owned = max((num(rows[s]["owned_bp"], int) or 0 for s in tc_samples if s in rows), default=0)
            for s in tc_samples:
                r = rows.get(s)
                if not r:
                    continue
                owned, frac = num(r["owned_bp"], int) or 0, num(r["frac_in_trc"]) or 0.0
                # >= min_trc_cov of its bp inside TideCluster arrays is array
                # evidence by itself, tandem_family flag or not (long units and
                # few copies per array can miss family_tandem's thresholds)
                ok = frac >= t.min_trc_cov
                # A tandem family whose repeat unit is longer than TideHunter looked
                # for can't be in a TRC: no TideCluster evidence either way.
                ftr = ft.get(fam, {}).get(s, {})
                unit = kite_unit or num(ftr.get("monomer_period")) or num(ftr.get("cons_len")) or 0
                too_long = (not ok and r["tandem_family"] == "True"
                            and (num(ftr.get("tandem_frac")) or 0) >= t.min_tandem_frac_long
                            and unit > max_period.get(s, 3000))
                share_ok = r["tandem_family"] == "True" or owned >= t.eligible_min_share * top_owned
                if owned >= t.major_min_bp * scale.get(s, 1) and share_ok:
                    elig.append(s)
                    if ok:
                        passed.append(s)
                    elif too_long:
                        beyond.append(s)
                    elif r["tandem_family"] == "True" and frac >= t.partial_trc_cov:
                        partial.append(s)
                elif ok or frac >= t.partial_trc_cov:
                    weak.append(s)
            inferred = []
            for target, donor in infer.items():
                r = rows.get(target)
                if not r:
                    continue
                if donor in passed:
                    has = r["tandem_family"] == "True" and (num(r["owned_bp"], int) or 0) >= t.major_min_bp * scale.get(target, 1)
                    inferred.append(f"{target}:{'confirms' if has else 'no_arrays'}")
            b_trc = best(cov_lib, fam, "trc")
            b_rd = best(cov_rdna, fam, "rdna")
            b_mt = best(cov_mito, fam, "mito")
            rd_flag = any(r["trc_rdna_flag"] == "True" for r in rows.values())

            prop, conf, why = "NA", "none", []
            if cf.split("/")[0] in ("Simple_repeat", "Low_complexity"):
                why.append("simple repeat (RepeatMasker built-in): not called")
            elif b_rd and b_rd[0] >= t.rdna_min_cov:
                prop, conf = "rRNA", "high"
                why.append(f"{b_rd[0]:.2f} of consensus matches ribotin {strip(b_rd[2])} at {b_rd[1]:.1f}%")
            elif b_mt and b_mt[0] >= t.mito_min_cov:
                prop, conf = "Other/NUMT", "high"
                why.append(f"{b_mt[0]:.2f} of consensus matches mitogenome {strip(b_mt[2])} at {b_mt[1]:.1f}%")
            elif elig and passed and len(passed) == len(elig):
                prop, conf = "Satellite", "high"
                why.append(f"tandem and >= {t.min_trc_cov} of its bp in TRC arrays in all eligible "
                           f"samples ({','.join(passed)})")
            elif elig and beyond and len(passed) + len(beyond) == len(elig):
                # TE-labelled families in arrays may be tandem segmental copies of
                # the TE: RepeatMasker arrays alone only earn them "low".
                prop, conf = "Satellite", "medium" if cf in ("Unknown", "NA") or cf.startswith("Satellite") else "low"
                lim = ",".join(f"{s}:-P {max_period.get(s, 3000)}" for s in beyond)
                found = ",".join(
                    f"{s}:{(num(rows[s]['tidehunter_bp'], int) or 0) / max(num(rows[s]['owned_bp'], int) or 1, 1):.2f}"
                    for s in beyond)
                why.append(f"unit longer than TideHunter's detection limit ({lim}): its arrays are only found "
                           f"where a sub-period <= that exists (share of its bp in TideHunter arrays: {found}); "
                           f"tandem by RepeatMasker arrays (family_tandem)"
                           + (f"; TRC support in {','.join(passed)}" if passed else "")
                           + ("" if conf == "medium" else "; tandemly arrayed TE (segmental copies?)"))
            elif elig and partial and len(passed) + len(partial) + len(beyond) == len(elig):
                prop, conf = "Satellite", "medium"
                fr = ",".join(f"{s}:{rows[s]['frac_in_trc']}" for s in partial)
                why.append(f"partial TRC support (frac_in_trc {fr}; >= {t.partial_trc_cov}, < {t.min_trc_cov})"
                           + (f"; full TRC support in {','.join(passed)}" if passed else ""))
            elif passed or weak:
                prop, conf = "Satellite", "medium" if passed else "low"
                why.append(f"TRC support in {','.join(passed + weak)} of eligible {','.join(elig) or 'none'}")
            if prop in ("Satellite",) and rd_flag:
                prop, conf = "rRNA", "medium"
                why.append("its top TRC is flagged rDNA by TideCluster (no ribotin match)")
            tier, m = motif_tier(fam)
            if tier and prop not in ("rRNA", "Other/NUMT") and cf.split("/")[0] not in ("Simple_repeat", "Low_complexity"):
                c_, i_, name = m
                rank = {"none": 0, "low": 1, "medium": 2, "high": 3}
                if tier == "partial":
                    why.append(f"partial match to harmonized motif {name} ({c_:.2f} of consensus at {i_:.1f}%): "
                               f"a satellite segment plus other sequence (chimeric consensus?)")
                    if prop == "NA":
                        prop, conf = "Satellite", "low"
                else:
                    why.append(f"{c_:.2f} of consensus matches harmonized motif {name} at {i_:.1f}%")
                    # The harmonized library holds TE-derived "satellite" units
                    # too, so a motif match alone does not make a dispersed TE a
                    # satellite: it needs array evidence in some sample (unless
                    # already labelled Satellite), and a TE-labelled family stays
                    # at most medium (TE-derived satellite).
                    arrays = [s for s, r in rows.items()
                              if r["tandem_family"] == "True"
                              or (num(ft.get(fam, {}).get(s, {}).get("tandem_frac")) or 0) >= t.motif_array_tandem_frac
                              or (num(r["frac_in_trc"]) or 0) >= t.partial_trc_cov]
                    te = cf.split("/")[0] in TE_CLASSES
                    sat_label = cf.startswith("Satellite")   # labelled already: the motif only names it
                    cap = tier if sat_label or (arrays and not te) else "medium" if arrays else "low"
                    if not arrays and not sat_label:
                        why.append("no array evidence (not tandem, not in TRC arrays): TE-derived motif? capped at low")
                    elif te and tier == "high":
                        why.append(f"labelled {cf}: TE-derived satellite? capped at medium")
                    tier_c = min(tier, cap, key=rank.get)
                    if prop in ("NA", "Satellite") and rank[tier_c] > rank[conf]:
                        prop, conf = "Satellite", tier_c
            if prop != "NA" and cf not in ("NA", "Unknown", prop) and not cf.startswith(prop):
                why.append(f"currently {cf} ({cf}-derived)")
            calls[fam] = {"class_family": cf, "prop": prop, "conf": conf, "motif": m if tier else None,
                          "motif_tier": tier}
            out.write(f"{fam}\t{cf}\t{len(passed)}\t{len(elig)}\t{','.join(beyond) or '.'}\t{','.join(inferred) or '.'}\t"
                      + (f"{strip(b_trc[2])}:cov={b_trc[0]:.2f}:id={b_trc[1]:.1f}" if b_trc else "NA") + "\t"
                      + (f"{m[2]}\t{m[0]:.2f}\t{m[1]:.1f}\t{tier or 'below'}" if m else "NA\tNA\tNA\tNA") + "\t"
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
        """Repeat unit: KITE founder (sample with the most TRC bp), else the
        TideHunter median (sample with the most TideHunter bp), else the
        consensus self-period."""
        def top(col):
            return max(tc_samples, key=lambda s: num(xc.get(fam, {}).get(s, {}).get(col), int) or 0, default=None)
        s_k = top("in_trc_bp")
        m = num(xc.get(fam, {}).get(s_k, {}).get("kite_founder_median")) if s_k else None
        if m:
            return m, f"kite:{s_k}"
        best_s = top("tidehunter_bp")
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

    # ---- TRC consensus matches across samples (rotation-safe: TRC consensi are dimers)
    def trc_name(x):
        sample, trc = strip(x).split(":", 1)
        return sample, trc

    trc_match = {}
    trc_rows = []
    seen = set()
    for (q, sub), (c, i) in cov_lib.items():
        if not (q.startswith("trc" + SEP) and sub.startswith("trc" + SEP)):
            continue
        (sq, tq), (ss_, ts) = trc_name(q), trc_name(sub)
        if sq == ss_:
            continue
        key = tuple(sorted(((sq, tq), (ss_, ts))))
        if key in seen:
            continue
        seen.add(key)
        (s1, t1), (s2, t2) = key
        c12 = cov_lib.get((f"trc{SEP}{s1}:{t1}", f"trc{SEP}{s2}:{t2}"), (0.0, 0.0))
        c21 = cov_lib.get((f"trc{SEP}{s2}:{t2}", f"trc{SEP}{s1}:{t1}"), (0.0, 0.0))
        ident = max(c12[1], c21[1])
        if max(c12[0], c21[0]) < t.min_pair_cov or ident < t.related_min_id:
            continue
        f1 = num(trc_info.get((s1, t1), {}).get("kite_founder_median"))
        f2 = num(trc_info.get((s2, t2), {}).get("kite_founder_median"))
        mm = monomer_match(f1, f2, t.monomer_tol)
        trc_match[key] = (max(c12[0], c21[0]), ident, mm)
        trc_rows.append((s1, t1, s2, t2, c12[0], c21[0], ident, f1, f2, mm))
    with open(t.trc_pairs_out, "w") as out:
        out.write("sample_a\ttrc_a\tsample_b\ttrc_b\tcov_a_by_b\tcov_b_by_a\tidentity\t"
                  "kite_founder_a\tkite_founder_b\tmonomer_match\n")
        for s1, t1, s2, t2, a_, b_, i_, f1, f2, mm in sorted(trc_rows):
            out.write(f"{s1}\t{t1}\t{s2}\t{t2}\t{a_:.2f}\t{b_:.2f}\t{i_:.1f}\t"
                      f"{'NA' if f1 is None else f'{f1:g}'}\t{'NA' if f2 is None else f'{f2:g}'}\t{mm}\n")

    def top_trcs(fam):
        return {s: xc[fam][s]["top_trc"] for s in tc_samples
                if xc.get(fam, {}).get(s, {}).get("top_trc") not in (None, "NA")}

    def top_trc_match(a, b):
        """Same top TRC in a sample, or matching top-TRC consensi across samples."""
        ta, tb = top_trcs(a), top_trcs(b)
        hits = []
        for sa_, x in ta.items():
            for sb_, y in tb.items():
                if sa_ == sb_ and x == y:
                    hits.append(f"{sa_}:{x}=same")
                elif sa_ != sb_:
                    m = trc_match.get(tuple(sorted(((sa_, x), (sb_, y)))))
                    if m:
                        hits.append(f"{sa_}:{x}~{sb_}:{y}(cov={m[0]:.2f},id={m[1]:.1f},{m[2]})")
        return hits

    def kite_in(fam, s):
        return num(xc.get(fam, {}).get(s, {}).get("kite_founder_median"))

    with open(t.pairs_out, "w") as out:
        out.write("family_a\tfamily_b\tcov_a_by_b\tcov_b_by_a\tidentity\tshared_trc_frac\t"
                  "max_shared_trc_frac\tmonomer_a\tmonomer_b\tmonomer_match\tsame_trc_superfamily\t"
                  "top_trc_match\tverdict\tnote\n")
        for a, b in sorted(pairs):
            ab = cov_lib.get((f"lib{SEP}{a}", f"lib{SEP}{b}"))
            ba = cov_lib.get((f"lib{SEP}{b}", f"lib{SEP}{a}"))
            cov_ab, cov_ba = (ab[0] if ab else 0.0), (ba[0] if ba else 0.0)
            ident = max([x[1] for x in (ab, ba) if x], default=0.0)
            shared, per, by_s = [], [], {}
            for s in tc_samples:
                ta, tb = trc_bp(a, s), trc_bp(b, s)
                if not ta or not tb:
                    per.append(f"{s}:NA")
                    continue
                # A family with only a sprinkle of bp in this sample's TRCs would
                # score 1.00 (overlap coefficient); require real arrays in both.
                if min(sum(ta.values()), sum(tb.values())) < t.min_shared_trc_bp:
                    per.append(f"{s}:low_bp({sum(ta.values())}/{sum(tb.values())})")
                    continue
                v = sum(min(ta[k], tb[k]) for k in set(ta) & set(tb)) / min(sum(ta.values()), sum(tb.values()))
                shared.append(v)
                by_s[s] = v
                per.append(f"{s}:{v:.2f}")
            mx = max(shared, default=None)
            # Monomers from the sample where they share arrays most (both need a
            # KITE founder there); else each family's own best estimate.
            both = [s for s in sorted(by_s, key=lambda x: -by_s[x]) if kite_in(a, s) and kite_in(b, s)]
            if both:
                ma, mb = kite_in(a, both[0]), kite_in(b, both[0])
                sa = sb = f"kite:{both[0]}"
            else:
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
            elif co and mono_ok and mm != "NA":
                # consensi don't align; matching TideCluster consensi across samples
                # make it the same satellite, otherwise check by hand
                cross = [h for h in top_trc_match(a, b)
                         if "~" in h and (",same)" in h or ",multiple_" in h)]
                verdict = "same_satellite_diverged" if cross else "same_arrays_same_monomer"
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
                      f"{mm}\t{same_sf}\t{';'.join(top_trc_match(a, b)) or '.'}\t{verdict}\t"
                      f"{'; '.join(note) or '.'}\n")

    # ---- proposals (curated_families format)
    with open(t.proposals_out, "w") as out:
        out.write("group\tfamily\tclass_family\tpart\tevidence\tnote\tconfidence\n")
        for fam in families:
            c = calls[fam]
            if c["conf"] not in ("high", "medium"):
                continue
            # Named after the matching harmonized motif (this family's, a partial
            # match included, else its same_satellite group's high/medium one);
            # otherwise same_satellite pairs share one group name (never
            # merged, only grouped)
            root = find(fam)
            named = next((calls[x]["motif"][2] for x in (fam, root)
                          if x in calls and calls[x]["motif"]
                          and calls[x]["motif_tier"] in (("high", "medium", "partial") if x == fam
                                                         else ("high", "medium"))), None)
            group = named if named and c["prop"] == "Satellite" else f"{c['prop'].split('/')[-1]}-{root}"
            note = f"{c['class_family']}-derived" if c["class_family"] not in ("NA", "Unknown", c["prop"]) else "."
            part = "45S_unit" if c["prop"] == "rRNA" else "partial" if c["motif_tier"] == "partial" else "."
            out.write(f"{group}\t{fam}\t{c['prop']}\t{part}\tsatellite_crosscheck\t{note}\t{c['conf']}\n")


def cmd_apply(args):
    from family_groups import read_groups
    user = {fam for _g, fam, _r in read_groups(args.user_table)} if args.user_table else set()
    keep = set(args.confidence)
    n_in = n_out = n_user = 0
    with open(args.proposals) as fh, open(args.out, "w") as out:
        header = fh.readline()
        out.write(header)
        cols = header.rstrip("\n").split("\t")
        i_fam, i_conf = cols.index("family"), cols.index("confidence")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            n_in += 1
            if f[i_conf] not in keep:
                continue
            if f[i_fam] in user:
                n_user += 1
                continue
            out.write(line)
            n_out += 1
    print(f"[satellite_evidence] apply: {n_out} of {n_in} proposals applied (confidence {sorted(keep)}); "
          f"{n_user} left to the user's curated table", file=sys.stderr)


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
    s.add_argument("--motifs", default="", help="harmonized satellite motif library (FASTA, one unit each)")
    s.add_argument("--motif-tile", type=int, default=300, help="tile each motif to at least this many bp")
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
    r.add_argument("--partial-trc-cov", type=float, default=0.3)
    r.add_argument("--min-shared-trc-bp", type=int, default=50000)
    r.add_argument("--copy-scale", nargs="*", default=[], help="sample=N: multiply --major-min-bp (2 for dual_hap)")
    r.add_argument("--tc-params", nargs="*", default=[], help="tidecluster_regions.py --params tables")
    r.add_argument("--min-shared-trc-frac", type=float, default=0.5)
    r.add_argument("--max-independent-shared-frac", type=float, default=0.1)
    r.add_argument("--monomer-tol", type=float, default=0.05)
    r.add_argument("--rdna-min-cov", type=float, default=0.5)
    r.add_argument("--rdna-min-id", type=float, default=90.0)
    r.add_argument("--mito-min-cov", type=float, default=0.5)
    r.add_argument("--mito-min-id", type=float, default=80.0)
    r.add_argument("--eligible-min-share", type=float, default=0.1,
                   help="a non-tandem sample counts as eligible only with >= this share of the family's top-sample bp")
    r.add_argument("--motif-min-id", type=float, default=75.0, help="HSP identity floor for motif coverage")
    r.add_argument("--motif-high-cov", type=float, default=0.7)
    r.add_argument("--motif-high-id", type=float, default=85.0)
    r.add_argument("--motif-min-cov", type=float, default=0.5, help="medium: this coverage at >= motif-min-id")
    r.add_argument("--motif-array-tandem-frac", type=float, default=0.5,
                   help="a motif match counts above low only when some sample shows arrays: tandem_family, "
                        "tandem_frac >= this, or frac_in_trc >= --partial-trc-cov")
    r.add_argument("--motif-partial-cov", type=float, default=0.15,
                   help="partial (chimera flag): this coverage at >= motif-high-id")
    r.add_argument("--evidence-out", required=True)
    r.add_argument("--calls-out", required=True)
    r.add_argument("--pairs-out", required=True)
    r.add_argument("--proposals-out", required=True)
    r.add_argument("--trc-pairs-out", required=True)

    a = sub.add_parser("apply")
    a.add_argument("--proposals", required=True)
    a.add_argument("--user-table", default="", help="classify.curated_families (wins on conflict)")
    a.add_argument("--confidence", nargs="+", default=["high"])
    a.add_argument("--out", required=True)

    args = ap.parse_args()
    {"select": cmd_select, "report": cmd_report, "apply": cmd_apply}[args.cmd](args)


if __name__ == "__main__":
    main()
