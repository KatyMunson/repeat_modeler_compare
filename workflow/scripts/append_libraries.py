#!/usr/bin/env python3
"""Assemble a masking library: a base FASTA (clustered shared de novo
library, one species' prefixed families, or a curated override), then an
optional Dfam export appended as-is. Fails on a family name present in
both (a curated override that already contains Dfam entries should not
get them twice silently). Writes --report (source, path, n_appended).
Stdlib only."""

import argparse
import sys

from fasta_utils import iter_fasta, write_fasta


def name_of(header):
    return header.split()[0].split("#", 1)[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base", required=True)
    ap.add_argument("--dfam")
    ap.add_argument("--out", required=True)
    ap.add_argument("--report", required=True)
    args = ap.parse_args()

    report = []
    base_names = set()
    with open(args.out, "w") as out:
        n = 0
        for header, seq in iter_fasta(args.base):
            base_names.add(name_of(header))
            write_fasta(out, header, seq)
            n += 1
        report.append(("base", args.base, n))

        if args.dfam:
            n = 0
            for header, seq in iter_fasta(args.dfam):
                if name_of(header) in base_names:
                    raise ValueError(f"{args.dfam}: Dfam family '{name_of(header)}' is already in {args.base}")
                write_fasta(out, header, seq)
                n += 1
            report.append(("dfam", args.dfam, n))

    with open(args.report, "w") as fh:
        fh.write("source\tpath\tn_appended\n")
        for row in report:
            fh.write("\t".join(str(x) for x in row) + "\n")
    print("[append_libraries] " + "; ".join(f"{s}:{n}" for s, _, n in report), file=sys.stderr)


if __name__ == "__main__":
    main()
