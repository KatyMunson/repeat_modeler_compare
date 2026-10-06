#!/usr/bin/env python3
"""Manifest v2 parser, shared by the Snakefile and combine_summaries.py.

Header-driven: the first non-'#' line names the columns (any order). One row
per assembly (sample). Every cell must hold a value with no spaces; "not
known / not applicable" is written as one of MISSING (case-insensitive).

  sample_id      [A-Za-z0-9]+, unique; wildcard and library-name prefix
  taxon          [A-Za-z0-9_.-]+, e.g. Sturnella_magna; samples of one species
                 share it (cluster sharing per taxon, cross-check inference)
  taxid          NCBI taxid (digits) or missing
  fasta          path to the assembly FASTA (may be gzipped); must exist
  tissue         germline | soma | unknown
  sex            ZW | ZZ | XX | XY | unknown
  assembly_type  haploid | primary | hap1 | hap2 | dual_hap | unknown
                 (dual_hap: both haplotypes in one assembly, e.g. unphased
                 Verkko; absolute copy-number thresholds are doubled)
  accession      free token or missing

Stdlib only."""

import os
import re
import sys

COLUMNS = ("sample_id", "taxon", "taxid", "fasta", "tissue", "sex", "assembly_type", "accession")
MISSING = {".", "na", "no", "false", "none", "unknown"}
ENUMS = {
    "tissue": ("germline", "soma", "unknown"),
    "sex": ("ZW", "ZZ", "XX", "XY", "unknown"),
    "assembly_type": ("haploid", "primary", "hap1", "hap2", "dual_hap", "unknown"),
}
# Haplotype copies per locus, by assembly_type: scales absolute copy thresholds.
COPIES = {"dual_hap": 2}


def is_missing(value):
    return value.strip().lower() in MISSING


def _legacy_hint(path, rows):
    example = ["\t".join(COLUMNS)]
    for f in rows[:4]:
        if len(f) >= 5:
            sid, name, fasta, tissue, acc = f[:5]
            taxon = re.sub(r"[^A-Za-z0-9_.-]", "_", name.strip()) or "NA"
            example.append("\t".join([sid, taxon, "NA", fasta, tissue, "unknown", "unknown", acc or "NA"]))
    return (f"{path} is a legacy 5-column manifest (species_id species_name fasta tissue accession). "
            f"Manifest v2 needs a header line with the columns {', '.join(COLUMNS)} "
            f"(see manifest.tsv). Your rows converted (fill in taxid, sex, assembly_type):\n  "
            + "\n  ".join(example))


def parse_manifest(path, check_fasta=True):
    """List of dicts (one per sample, manifest order) with the COLUMNS keys,
    plus 'copies' (1, or 2 for dual_hap). Raises ValueError with line numbers."""
    header, rows, raw = None, [], []
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n").rstrip("\r")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split("\t")
            if header is None:
                if "sample_id" not in fields:
                    raw.append(fields)
                    for l2 in fh:
                        if l2.strip() and not l2.startswith("#"):
                            raw.append(l2.rstrip("\n").split("\t"))
                    raise ValueError(_legacy_hint(path, raw))
                missing = [c for c in COLUMNS if c not in fields]
                if missing:
                    raise ValueError(f"{path}:{lineno}: header lacks column(s) {missing}; "
                                     f"manifest v2 columns are {', '.join(COLUMNS)}")
                header = fields
                continue
            if len(fields) != len(header):
                raise ValueError(f"{path}:{lineno}: {len(fields)} fields, header has {len(header)}")
            row = dict(zip(header, fields))
            for col in COLUMNS:
                v = row[col]
                if v == "" or v != v.strip() or " " in v:
                    raise ValueError(f"{path}:{lineno}: column '{col}' is blank or contains spaces "
                                     f"(write NA / . / unknown for 'not known')")
            rows.append((lineno, row))
    if header is None:
        raise ValueError(f"{path}: no header line found")

    out, seen = [], set()
    for lineno, row in rows:
        sid = row["sample_id"]
        if not re.fullmatch(r"[A-Za-z0-9]+", sid):
            raise ValueError(f"{path}:{lineno}: sample_id '{sid}' must match [A-Za-z0-9]+")
        if sid in seen:
            raise ValueError(f"{path}:{lineno}: duplicate sample_id '{sid}'")
        seen.add(sid)
        if is_missing(row["taxon"]) or not re.fullmatch(r"[A-Za-z0-9_.-]+", row["taxon"]):
            raise ValueError(f"{path}:{lineno}: taxon '{row['taxon']}' must be a name like Sturnella_magna "
                             f"([A-Za-z0-9_.-]+)")
        if not is_missing(row["taxid"]) and not row["taxid"].isdigit():
            raise ValueError(f"{path}:{lineno}: taxid '{row['taxid']}' must be digits or NA")
        for col, allowed in ENUMS.items():
            v = row[col]
            canon = next((a for a in allowed if a.lower() == v.lower()), None)
            if canon is None and is_missing(v):
                canon = "unknown"
            if canon is None:
                raise ValueError(f"{path}:{lineno}: {col} must be one of {'|'.join(allowed)}, got '{v}'")
            row[col] = canon
        if check_fasta and not os.path.exists(row["fasta"]):
            raise ValueError(f"{path}:{lineno}: fasta path does not exist: {row['fasta']}")
        rec = {c: row[c] for c in COLUMNS}
        rec["copies"] = COPIES.get(rec["assembly_type"], 1)
        out.append(rec)
    if len(out) < 2:
        raise ValueError(f"{path}: at least 2 samples are required, found {len(out)}")
    return out


if __name__ == "__main__":
    for r in parse_manifest(sys.argv[1], check_fasta=False):
        print(r)
