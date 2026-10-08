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
Stdlib only."""

import argparse
import csv
import os
import sys

from family_groups import curated_classes

EVAL_RANK = {"Perfect": 3, "Good": 2, "Reco_check": 1, "Need_check": 0, "Need_ext": 0}
EVAL_CONF = {"Perfect": "high", "Good": "medium"}


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

    picked = []
    for fam, by_s in owned.items():
        if fam in tandem or not unresolved(label.get(fam, "Unknown")):
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
          f"(top {args.top_n}, >= {args.min_bp:,} bp)", file=sys.stderr)


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
                rows.append(r)

    cols = ["sample", "family", "input_label", "output_name", "status", "evaluation", "low_copy",
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
            best = max(outs, key=lambda r: (EVAL_RANK.get(r["evaluation"], 0), num(r["output_genome_cov_len"])))
            new = best.get("output_TE_type", "")
            if not new or new.lower() in ("nan", "unknown") or new == best["input_label"]:
                continue
            conf = EVAL_CONF.get(best["evaluation"], "low")
            n_conf[conf] += 1
            note = [f"TEtrimmer {best['evaluation']} {best['output_name']} "
                    f"({best['input_length']} -> {best['output_length']} bp)"]
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
    args = ap.parse_args()
    {"select": cmd_select, "report": cmd_report}[args.cmd](args)


if __name__ == "__main__":
    main()
