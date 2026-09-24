#!/usr/bin/env python3
"""Assemble a masking library: a base FASTA (clustered shared de novo
library, one species' prefixed families, or a curated override), then an
optional Dfam export, then zero or more satellite libraries -- appended
as-is, never clustered or prefixed (their names are already
taxonomy-coded by compare_assemblies_satellites stage 02b).

Satellite library checks (fail loudly, naming the offending entry):
  - every header must be "<name>#<Class/Family>"
  - the same library passed twice (by content md5) is appended once
  - a satellite name that already exists in the base library is skipped
    (e.g. a curated_override that already contains it) and reported
  - the same name with a different sequence across libraries, or a
    satellite name colliding with a Dfam name, is an error
Writes --report (source, path, n_appended, n_skipped). Stdlib only."""

import argparse
import hashlib
import re
import sys

from fasta_utils import iter_fasta, write_fasta

HEADER = re.compile(r"^[^#\s]+#[^\s#]+$")


def name_of(header):
    return header.split()[0].split("#", 1)[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base", required=True)
    ap.add_argument("--dfam")
    ap.add_argument("--satellite-libs", nargs="*", default=[])
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
        report.append(("base", args.base, n, 0))

        dfam_names = set()
        if args.dfam:
            n = 0
            for header, seq in iter_fasta(args.dfam):
                dfam_names.add(name_of(header))
                write_fasta(out, header, seq)
                n += 1
            report.append(("dfam", args.dfam, n, 0))

        seen_md5 = set()
        sat_seq = {}
        for path in args.satellite_libs:
            with open(path, "rb") as fh:
                md5 = hashlib.md5(fh.read()).hexdigest()
            if md5 in seen_md5:
                report.append(("satellite_duplicate_file", path, 0, 0))
                continue
            seen_md5.add(md5)
            n = skipped = 0
            for header, seq in iter_fasta(path):
                token = header.split()[0]
                if not HEADER.match(token):
                    raise ValueError(f"{path}: satellite header '{header}' is not '<name>#<Class/Family>'")
                name = name_of(token)
                if name in dfam_names:
                    raise ValueError(f"{path}: satellite name '{name}' collides with a Dfam family name")
                if name in sat_seq:
                    if sat_seq[name] != seq.upper():
                        raise ValueError(f"{path}: satellite name '{name}' already appended with a different sequence")
                    skipped += 1
                    continue
                if name in base_names:
                    skipped += 1
                    continue
                sat_seq[name] = seq.upper()
                write_fasta(out, token, seq)
                n += 1
            report.append(("satellite", path, n, skipped))

    with open(args.report, "w") as fh:
        fh.write("source\tpath\tn_appended\tn_skipped\n")
        for row in report:
            fh.write("\t".join(str(x) for x in row) + "\n")
    print("[append_libraries] " + "; ".join(f"{s}:{a}+{k}skip" for s, _, a, k in report), file=sys.stderr)


if __name__ == "__main__":
    main()
