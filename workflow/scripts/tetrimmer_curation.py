#!/usr/bin/env python3
"""Feed the largest unresolved Unknown families to TEtrimmer and read its
results back as curated-table proposals (README "TEtrimmer").

Subcommands:
  select   per sample, the families still Unknown (or <Order>/Unknown)
           after the curated, cross-check, reclassification and
           family_neighbors labels, not tandem in any sample, whose largest
           footprint (owned_bp) is in this sample: the top --top-n by
           owned_bp with >= --min-bp. Each family runs once, against the
           genome where it is most abundant. Writes their consensi as
           >family#class for TEtrimmer --input_file.
  report   read each sample's TEtrimmer output directory
           (Sequence_name_mapping.txt, summary.txt), map TEtrimmer's
           sanitized names back to library family IDs, and write
             --summary-out    one row per TEtrimmer output consensus
             --proposals-out  classify.curated_families format, one row per
                              family whose label TEtrimmer changed: its best
                              output consensus (evaluation, then genome
                              coverage); confidence Perfect -> high,
                              Good -> medium, else low
           Report-only: review the PDFs under TEtrimmer_for_proof_curation/
           and copy accepted rows into your curated table.
  consensi TEtrimmer's final consensi (TEtrimmer_consensus.fasta) for every
           processed output, one FASTA for the project + an id map.
  tesorter TEsorter calls on those extended consensi against the calls on
           the original library consensus (classify step) and the
           classifier label (--classified: RepeatClassifier in the TEtools
           container; else TEtrimmer's own): --out one row per output
           consensus (domain gained / kept / changed / lost / none; agrees
           with the classifier or not); --proposals-out one row per
           still-unresolved family. A domain call is proposed first: high =
           TEtrimmer Perfect/Good and the classifier gives the same
           superfamily, low = the classifier says another order, else
           medium. Without a domain, a resolved classifier label: medium
           for Perfect/Good, else low. A call of a different order than
           the family's existing <Order>/Unknown label (e.g. LINE for an
           LTR-pipeline LTR/Unknown family) is always low, noted
           "ORDER CHANGE".
Stdlib only."""

import argparse
import csv
import os
import sys

from family_groups import curated_classes

EVAL_RANK = {"Perfect": 3, "Good": 2, "Reco_check": 1, "Need_check": 0, "Need_ext": 0}
EVAL_CONF = {"Perfect": "high", "Good": "medium"}
# TEtrimmer 1.7.2 means to promote Reco_check to Good at >= 3 full-length
# genomic copies (boundarycrop.py), but tests for "Recon_check", so it never
# does. The reports apply that rule as evaluation_regraded.
PROMOTE_MIN_FULL = 3


def regrade(evaluation, full_n):
    """(evaluation used by the reports, note): Reco_check with >=
    PROMOTE_MIN_FULL full-length copies -> Good."""
    try:
        n = int(float(full_n))
    except (TypeError, ValueError):
        n = 0
    if evaluation == "Reco_check" and n >= PROMOTE_MIN_FULL:
        return "Good", f"regraded Reco_check -> Good ({n} full-length copies; TEtrimmer promotion bug)"
    return evaluation, ""


def bare(name):
    return name.split("#", 1)[0]


def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            yield dict(zip(header, line.rstrip("\n").split("\t")))


def unresolved(cls):
    return cls == "Unknown" or cls.endswith("/Unknown")


def cmd_select(args):
    owned, tandem, label = {}, set(), {}
    for r in read_tsv(args.family_tandem):
        fam = bare(r["family"])
        owned.setdefault(fam, {})[r["sample"]] = int(r["owned_bp"])
        label.setdefault(fam, r["class_family"])
        if r.get("tandem_family") == "True":
            tandem.add(fam)
    if args.reclass:
        for r in read_tsv(args.reclass):
            if r.get("new_class") and r["new_class"] != r.get("old_class"):
                label[bare(r["family"])] = r["new_class"]
    for path in args.curated:
        label.update(curated_classes(path))
    if args.neighbor_proposals:
        for r in read_tsv(args.neighbor_proposals):
            if r.get("confidence") in args.neighbor_confidence:
                label[bare(r["family"])] = r["class_family"]

    # Only families with a sequence in --library (de novo consensi; Dfam
    # entries are left out of it unless classify.screen_dfam), so a curated
    # Dfam Unknown such as UCON1 cannot take a top_n slot.
    with open(args.library) as fh:
        in_lib = {bare(l[1:].split()[0]) for l in fh if l.startswith(">")}
    picked, n_other = [], 0
    for fam, by_s in owned.items():
        if fam in tandem or not unresolved(label.get(fam, "Unknown")):
            continue
        if fam not in in_lib:
            n_other += 1
            continue
        lead = max(by_s, key=by_s.get)
        if lead == args.sample and by_s[lead] >= args.min_bp:
            picked.append((by_s[lead], fam))
    picked.sort(reverse=True)
    picked = picked[:args.top_n]
    want = {fam for _bp, fam in picked}

    n = 0
    with open(args.library) as fh, open(args.out, "w") as out:
        keep = False
        for line in fh:
            if line.startswith(">"):
                fam = bare(line[1:].split()[0])
                keep = fam in want
                if keep:
                    out.write(f">{fam}#{label.get(fam, 'Unknown')}\n")
                    want.discard(fam)
                    n += 1
            elif keep:
                out.write(line)
    if want:
        print(f"[tetrimmer select] WARNING {len(want)} families not in {args.library}: "
              f"{','.join(sorted(want)[:5])}...", file=sys.stderr)
    tot = sum(bp for bp, _f in picked)
    print(f"[tetrimmer select] {args.sample}: {n} families, {tot / 1e6:.1f} Mb owned "
          f"(top {args.top_n}, >= {args.min_bp:,} bp; {n_other} unresolved families not in "
          f"{os.path.basename(args.library)}, e.g. Dfam entries, not considered)", file=sys.stderr)


def cmd_report(args):
    rows = []
    for spec in args.run:
        sample, _, d = spec.partition("=")
        names = {}
        with open(os.path.join(d, "Sequence_name_mapping.txt")) as fh:
            fh.readline()
            for line in fh:
                orig, _, san = line.rstrip("\n").partition("\t")
                names[san] = orig
        with open(os.path.join(d, "summary.txt")) as fh:
            for r in csv.DictReader(fh):
                orig = names.get(r["input_name"], r["input_name"])
                r["family"], r["input_label"] = bare(orig), orig.partition("#")[2] or "Unknown"
                r["sample"] = sample
                r["evaluation_regraded"] = regrade(r.get("evaluation"), r.get("output_full_blast_n"))[0]
                rows.append(r)

    cols = ["sample", "family", "input_label", "output_name", "status", "evaluation", "evaluation_regraded",
            "low_copy",
            "input_length", "output_length", "in_out_identify", "output_TE_type", "output_terminal_repeat",
            "TSD", "start_pattern", "end_pattern", "database_hit", "input_genome_cov_len",
            "output_genome_cov_len", "output_full_blast_n", "cluster"]
    with open(args.summary_out, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join(str(r.get(c, "NA")).replace("\t", " ") for c in cols) + "\n")

    def num(x):
        try:
            return float(x)
        except (TypeError, ValueError):
            return 0.0

    by_fam = {}
    for r in rows:
        if r.get("status") != "processed" or r.get("evaluation") in ("NaN", "", None):
            continue
        by_fam.setdefault(r["family"], []).append(r)
    n_conf = {"high": 0, "medium": 0, "low": 0}
    with open(args.proposals_out, "w") as out:
        out.write("group\tfamily\tclass_family\tpart\tevidence\tnote\tconfidence\n")
        for fam in sorted(by_fam):
            outs = by_fam[fam]
            best = max(outs, key=lambda r: (EVAL_RANK.get(r["evaluation_regraded"], 0), num(r["output_genome_cov_len"])))
            new = best.get("output_TE_type", "")
            if not new or new.lower() in ("nan", "unknown") or new == best["input_label"]:
                continue
            conf = EVAL_CONF.get(best["evaluation_regraded"], "low")
            n_conf[conf] += 1
            note = [f"TEtrimmer {best['evaluation_regraded']} {best['output_name']} "
                    f"({best['input_length']} -> {best['output_length']} bp)"]
            rn = regrade(best.get("evaluation"), best.get("output_full_blast_n"))[1]
            if rn:
                note.append(rn)
            if best.get("output_terminal_repeat") not in ("False", "FALSE", "NaN", "", None):
                note.append(f"terminal repeat {best['output_terminal_repeat']}")
            if best.get("TSD") not in ("nan", "NaN", "", None, "False"):
                note.append(f"TSD {best['TSD']}")
            if best.get("database_hit") not in ("nan", "NaN", "", None, "False"):
                note.append(f"hits {best['database_hit']}")
            if len(outs) > 1:
                note.append(f"split into {len(outs)} consensi (subfamilies or chimera)")
            out.write(f"{fam}\t{fam}\t{new}\t.\ttetrimmer\t{'; '.join(note).replace(chr(9), ' ')}\t{conf}\n")
    print(f"[tetrimmer report] {len(rows)} output rows, {len(by_fam)} families processed; label changes: "
          f"{n_conf['high']} high, {n_conf['medium']} medium, {n_conf['low']} low", file=sys.stderr)


# ---------------------------------------------------------------- TEsorter check

def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            elif line:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def cmd_consensi(args):
    """TEtrimmer's final consensi (TEtrimmer_consensus.fasta, headers
    output_name#TE_type) for every processed row of tetrimmer_summary.tsv, as
    one FASTA with ids <sample>__<output_name>, plus a map id -> sample,
    family, output_name, TEtrimmer's final type."""
    rows = {}
    for r in read_tsv(args.summary):
        if r.get("status") == "processed" and r.get("output_name") not in ("", "NA", "NaN", None):
            rows[(r["sample"], r["output_name"])] = r
    runs = dict(spec.partition("=")[::2] for spec in args.run)
    n = missing = 0
    with open(args.out, "w") as out, open(args.map_out, "w") as mp:
        mp.write("id\tsample\tfamily\toutput_name\ttetrimmer_type\tlength\n")
        for sample, d in runs.items():
            path = os.path.join(d, "TEtrimmer_consensus.fasta")
            seqs = {}
            if os.path.exists(path):
                for h, sq in read_fasta(path):
                    name, _, te = h.partition("#")
                    seqs[name] = (te or "Unknown", sq)
            for (s, oname), r in rows.items():
                if s != sample:
                    continue
                if oname not in seqs:
                    missing += 1
                    continue
                te, sq = seqs[oname]
                cid = f"{sample}__{oname}"
                out.write(f">{cid}\n")
                for i in range(0, len(sq), 80):
                    out.write(sq[i:i + 80] + "\n")
                mp.write(f"{cid}\t{sample}\t{r['family']}\t{oname}\t{te}\t{len(sq)}\n")
                n += 1
    print(f"[tetrimmer consensi] {n} output consensi; {missing} processed outputs not found in "
          f"TEtrimmer_consensus.fasta", file=sys.stderr)


def _order(label):
    return (label or "").split("/", 1)[0]


def cmd_tesorter(args):
    """TEsorter on the extended consensi vs TEsorter on the original library
    consensus (classify step) vs TEtrimmer's own label. Table per output
    consensus; proposals (curated-table format) per family that is still
    unresolved and whose extended consensus carries a domain call."""
    from reclassify_unknown import domain_evidence, read_tesorter_cls
    from dna_te_pilot import same_superfamily

    rex, gydb = read_tesorter_cls(args.rexdb), read_tesorter_cls(args.gydb)
    lrex, lgydb = read_tesorter_cls(args.lib_rexdb), read_tesorter_cls(args.lib_gydb)
    summ = {(r["sample"], r["output_name"]): r for r in read_tsv(args.summary)}
    rc = {}
    if args.classified and os.path.exists(args.classified):
        for h, _sq in read_fasta(args.classified):
            cid, _, lab = h.partition("#")
            rc[cid] = lab or "Unknown"

    rows = []
    for m in read_tsv(args.map):
        s = summ.get((m["sample"], m["output_name"]), {})
        after, a_conf = domain_evidence(m["id"], rex, gydb, args.min_domains)
        before, b_conf = domain_evidence(m["family"], lrex, lgydb, args.min_domains)
        # the label to compare with: RepeatClassifier (TEtools) on the
        # extended consensus when given, else TEtrimmer's own
        tet = rc.get(m["id"], m["tetrimmer_type"]) if args.classified else m["tetrimmer_type"]
        if after and before:
            change = "kept" if _order(after) == _order(before) else "changed"
        elif after:
            change = "gained"
        elif before:
            change = "lost"
        else:
            change = "none"
        if not after:
            vs_tet = "."
        elif _order(tet) in ("Unknown", ""):
            vs_tet = "classifier_unknown"
        elif same_superfamily(after, tet) or after == tet:
            vs_tet = "agree"
        elif _order(after) == _order(tet):
            vs_tet = "agree_order"
        else:
            vs_tet = "disagree"
        rows.append({
            "sample": m["sample"], "family": m["family"], "input_label": s.get("input_label", "Unknown"),
            "output_name": m["output_name"], "evaluation": s.get("evaluation", "NA"),
            "evaluation_regraded": s.get("evaluation_regraded")
            or regrade(s.get("evaluation"), s.get("output_full_blast_n"))[0],
            "regrade_note": regrade(s.get("evaluation"), s.get("output_full_blast_n"))[1],
            "input_length": s.get("input_length", "NA"), "output_length": m["length"],
            "genome_cov": s.get("output_genome_cov_len", "NA"), "tetrimmer_type": m["tetrimmer_type"],
            "classifier_type": tet,
            "tesorter_original": before or ".", "tesorter_extended": after or ".",
            "domain_change": change, "vs_classifier": vs_tet,
            "conflict": "; ".join(x for x in (b_conf and f"original: {b_conf}", a_conf and f"extended: {a_conf}") if x) or ".",
        })

    cols = ["sample", "family", "input_label", "output_name", "evaluation", "evaluation_regraded",
            "input_length", "output_length",
            "genome_cov", "tetrimmer_type", "classifier_type", "tesorter_original", "tesorter_extended",
            "domain_change", "vs_classifier", "conflict"]
    with open(args.out, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join(str(r[c]).replace("\t", " ") for c in cols) + "\n")

    def num(x):
        try:
            return float(x)
        except (TypeError, ValueError):
            return 0.0

    # per family: its best output (evaluation, then coverage) carrying a
    # domain call, else a resolved classifier label
    def proposal(r):
        lab = r["input_label"]
        if not unresolved(lab):
            return None
        for call in (r["tesorter_extended"], r["classifier_type"]):
            if call in (".", "", None) or unresolved(call):
                continue
            if call == lab or ("/" not in call and _order(call) == _order(lab)):
                continue  # nothing beyond the label it already has (DNA for DNA/Unknown)
            return call
        return None

    best = {}
    for r in rows:
        call = proposal(r)
        if call is None:
            continue
        key = (r["tesorter_extended"] != ".", EVAL_RANK.get(r["evaluation_regraded"], 0), num(r["genome_cov"]))
        if r["family"] not in best or key > best[r["family"]][0]:
            best[r["family"]] = (key, r, call)
    n_conf = {"high": 0, "medium": 0, "low": 0}
    with open(args.proposals_out, "w") as out:
        out.write("group\tfamily\tclass_family\tpart\tevidence\tnote\tconfidence\n")
        for fam in sorted(best):
            _k, r, call = best[fam]
            good = r["evaluation_regraded"] in ("Perfect", "Good")
            domain = r["tesorter_extended"] != "."
            if not domain:
                conf = "medium" if good else "low"
                evidence = "tetrimmer+repeatclassifier"
            elif r["vs_classifier"] == "disagree":
                conf = "low"
                evidence = "tetrimmer+tesorter"
            elif r["vs_classifier"] == "agree" and good:
                conf = "high"
                evidence = "tetrimmer+tesorter+repeatclassifier"
            else:
                conf = "medium"
                evidence = "tetrimmer+tesorter"
            # the family already had an order (LTR/Unknown from the LTR pipeline,
            # DNA/Unknown, ...) and the call names another one: the earlier
            # label was wrong, or the extension drifted onto a nested /
            # neighbouring element. Review, never accept on the call alone.
            lab = r["input_label"]
            order_change = _order(lab) not in ("Unknown", "") and _order(call) != _order(lab)
            if order_change:
                conf = "low"
            n_conf[conf] += 1
            note = [f"ORDER CHANGE: contradicts the existing {lab}; check the PDF (earlier label wrong, "
                    f"or the extension drifted onto a nested / neighbouring element)" if order_change else "",
                    f"TEtrimmer {r['evaluation_regraded']} {r['output_name']} ({r['input_length']} -> {r['output_length']} bp)",
                    r["regrade_note"],
                    f"TEsorter {r['tesorter_extended']}" if domain else "no TEsorter domain",
                    f"RepeatClassifier {r['classifier_type']} ({r['vs_classifier']})" if domain
                    else f"RepeatClassifier {r['classifier_type']}",
                    {"gained": "no domain call on the original consensus",
                     "kept": f"original consensus also {r['tesorter_original']}",
                     "changed": f"original consensus called {r['tesorter_original']}"}.get(r["domain_change"], "")]
            if r["conflict"] != ".":
                note.append(r["conflict"])
            out.write(f"{fam}\t{fam}\t{call}\t.\t{evidence}\t"
                      f"{'; '.join(x for x in note if x).replace(chr(9), ' ')}\t{conf}\n")
    by = {}
    for r in rows:
        by[r["domain_change"]] = by.get(r["domain_change"], 0) + 1
    print(f"[tetrimmer tesorter] {len(rows)} extended consensi; domain calls {by}; proposals: "
          f"{n_conf['high']} high, {n_conf['medium']} medium, {n_conf['low']} low", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("select")
    s.add_argument("--family-tandem", required=True, help="combined family_tandem.tsv")
    s.add_argument("--library", required=True, help="library FASTA (family consensi)")
    s.add_argument("--sample", required=True)
    s.add_argument("--reclass", default="")
    s.add_argument("--curated", nargs="*", default=[], help="curated / applied tables")
    s.add_argument("--neighbor-proposals", default="", help="family_neighbors proposals")
    s.add_argument("--neighbor-confidence", nargs="*", default=["high"],
                   help="neighbor proposals at these confidences count as resolved")
    s.add_argument("--top-n", type=int, default=200)
    s.add_argument("--min-bp", type=int, default=500000)
    s.add_argument("--out", required=True)
    r = sub.add_parser("report")
    r.add_argument("--run", nargs="+", required=True, help="sample=TEtrimmer output directory")
    r.add_argument("--summary-out", required=True)
    r.add_argument("--proposals-out", required=True)
    c = sub.add_parser("consensi")
    c.add_argument("--run", nargs="+", required=True, help="sample=TEtrimmer output directory")
    c.add_argument("--summary", required=True, help="tetrimmer_summary.tsv (report)")
    c.add_argument("--out", required=True)
    c.add_argument("--map-out", required=True)
    t = sub.add_parser("tesorter")
    t.add_argument("--summary", required=True, help="tetrimmer_summary.tsv (report)")
    t.add_argument("--map", required=True, help="consensi --map-out")
    t.add_argument("--rexdb", required=True, help="TEsorter on the extended consensi")
    t.add_argument("--gydb", required=True)
    t.add_argument("--lib-rexdb", default="", help="TEsorter on the library consensi (classify step)")
    t.add_argument("--lib-gydb", default="")
    t.add_argument("--classified", default="", help="RepeatClassifier output on the extended consensi "
                   "(id#Class/Family); compared instead of TEtrimmer's own label")
    t.add_argument("--min-domains", type=int, default=1)
    t.add_argument("--out", required=True)
    t.add_argument("--proposals-out", required=True)
    args = ap.parse_args()
    {"select": cmd_select, "report": cmd_report, "consensi": cmd_consensi,
     "tesorter": cmd_tesorter}[args.cmd](args)


if __name__ == "__main__":
    main()
