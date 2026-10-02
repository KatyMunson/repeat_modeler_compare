#!/usr/bin/env python3
"""Build the TE-free host-protein set for the DIAMOND host-gene screen
(reclassify_unknown.py). Gene annotations still contain TE-derived models
(transposases, RT/Pol ORFs, "uncharacterized" ORFs of active TEs); left in,
they would relabel real TEs as host genes. Two filters:

  prep      keyword filter on FASTA descriptions (all inputs), IDs prefixed
            with the source file's stem so hits show the source species.
            Annotation proteomes go to --annot-out (to be domain-screened
            with TEsorter -st prot); Swiss-Prot goes to --sprot-out
            (keyword filter only: reviewed descriptions are reliable, and
            domain-screening ~570k proteins isn't worth it).
  finalize  drop annotation proteins TEsorter found TE domains in, then
            concatenate with the Swiss-Prot part.

Domesticated TE genes (PGBD, ZBED, CENP-B, ...) are removed too, on purpose:
a family matching one keeps its TE label from TEsorter instead of becoming
a "host gene". Counts are printed to stderr. Stdlib only."""

import argparse
import os
import sys

from fasta_utils import iter_fasta, write_fasta


def stem(path):
    base = os.path.basename(path)
    for ext in (".gz", ".faa", ".fasta", ".fa", ".pep", ".aa"):
        if base.endswith(ext):
            base = base[: -len(ext)]
    return base.replace(" ", "_")


def keyword_hit(desc, keywords):
    d = desc.lower()
    return any(k in d for k in keywords)


def filter_into(paths, keywords, out_fh, tag):
    for path in paths:
        kept = dropped = 0
        prefix = stem(path)
        for header, seq in iter_fasta(path):
            if keyword_hit(header, keywords):
                dropped += 1
                continue
            parts = header.split(None, 1)
            desc = parts[1] if len(parts) > 1 else ""
            write_fasta(out_fh, f"{prefix}__{parts[0]} {desc}".rstrip(), seq.rstrip("*"))
            kept += 1
        print(f"[host_proteins] {tag} {path}: kept {kept}, dropped {dropped} by TE keyword", file=sys.stderr)


def cmd_prep(args):
    keywords = [k.lower() for k in args.keyword]
    with open(args.annot_out, "w") as fh:
        filter_into(args.annotation, keywords, fh, "annotation")
    with open(args.sprot_out, "w") as fh:
        filter_into([args.swissprot] if args.swissprot else [], keywords, fh, "swissprot")


def te_domain_ids(paths):
    ids = set()
    for path in paths:
        if not os.path.exists(path):
            continue
        with open(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                tok = line.split("\t", 1)[0].split()[0]
                ids.add(tok)
                ids.add(tok.split("|", 1)[0])  # TEsorter may append |domain or |frame
    return ids


def cmd_finalize(args):
    te_ids = te_domain_ids(args.tesorter)
    kept = dropped = 0
    with open(args.out, "w") as out:
        for header, seq in iter_fasta(args.annot):
            pid = header.split()[0]
            if pid in te_ids:
                dropped += 1
                continue
            write_fasta(out, header, seq)
            kept += 1
        n_sprot = 0
        for header, seq in iter_fasta(args.sprot):
            write_fasta(out, header, seq)
            n_sprot += 1
    if te_ids and dropped == 0:
        print("[host_proteins] WARNING: TEsorter reported TE-domain proteins but none matched an "
              "annotation protein ID -- check ID formats", file=sys.stderr)
    print(f"[host_proteins] annotation proteins: kept {kept}, dropped {dropped} with TE domains; "
          f"Swiss-Prot proteins: {n_sprot}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prep")
    p.add_argument("--annotation", nargs="*", default=[], help="annotation proteomes (domain-screened next)")
    p.add_argument("--swissprot", default="")
    p.add_argument("--keyword", action="append", default=[])
    p.add_argument("--annot-out", required=True)
    p.add_argument("--sprot-out", required=True)
    f = sub.add_parser("finalize")
    f.add_argument("--annot", required=True)
    f.add_argument("--sprot", required=True)
    f.add_argument("--tesorter", nargs="*", default=[], help="TEsorter -st prot *.cls.tsv / *.dom.tsv")
    f.add_argument("--out", required=True)
    args = ap.parse_args()
    {"prep": cmd_prep, "finalize": cmd_finalize}[args.cmd](args)


if __name__ == "__main__":
    main()
