#!/usr/bin/env python3
"""Rewrite a species' RepeatModeler2 family FASTA headers from
'name#Class/Family [description]' to '{species_id}{sep}name#Class/Family',
preserving the '#Class/Family' suffix exactly and dropping any trailing
description. Fails loudly on duplicate names after prefixing. Stdlib only.

Optional protein-based host-gene filtering (config: library.protein_filter)
would plug in here, right after parsing and before writing the output FASTA
— e.g. BLASTing each family consensus against a protein database and
dropping hits above some identity/coverage threshold. Left as an explicit,
unimplemented hook (off by default) per the spec; not wired up in v1.
"""

import argparse
import sys


def iter_fasta(path):
    header = None
    seq_chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
    if header is not None:
        yield header, "".join(seq_chunks)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True, help="RepeatModeler2 *-families.fa for one species")
    ap.add_argument("--species-id", required=True)
    ap.add_argument("--sep", default="_")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    seen = set()
    n = 0
    with open(args.out, "w") as fh_out:
        for header, seq in iter_fasta(args.fasta):
            # Drop any trailing description after whitespace first.
            token = header.split()[0] if header.split() else header
            if "#" not in token:
                raise ValueError(
                    f"Family header '{header}' in {args.fasta} has no '#Class/Family' suffix"
                )
            name, class_family = token.split("#", 1)
            new_name = f"{args.species_id}{args.sep}{name}"

            if new_name in seen:
                raise ValueError(
                    f"Duplicate family name after prefixing: '{new_name}' "
                    f"(from original header '{header}' in {args.fasta})"
                )
            seen.add(new_name)

            fh_out.write(f">{new_name}#{class_family}\n")
            for i in range(0, len(seq), 60):
                fh_out.write(seq[i : i + 60] + "\n")
            n += 1

    print(f"[prefix_library] wrote {n} prefixed families to {args.out}", file=sys.stderr)
    if n == 0:
        raise ValueError(f"No families found in {args.fasta}")


if __name__ == "__main__":
    main()
