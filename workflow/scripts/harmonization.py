#!/usr/bin/env python3
"""Stage-02b (compare_assemblies_satellites) lookups used at DAG-build
time by the Snakefile. No NCBI calls: stage 02b's own cached lineages are
the source of truth for what its codes mean.

Species code resolution: taxonomy_cache/{code}.json holds
{"taxon": <the exact "Genus species" string stage 02b was given>, ...};
its file name is the collision-resolved species code 02b actually used
(e.g. EST, or EST2 after a collision). taxon_codes.tsv is deliberately NOT
used for species rows: it recomputes species codes lexically from the
NCBI scientific name, which can miss a stage-01 collision suffix or an
NCBI name change."""

import glob
import json
import os


def normalize_binomial(name):
    return " ".join(name.split()).lower()


def load_cache_codes(harmonization_dir):
    """{normalized binomial: code} from taxonomy_cache/*.json."""
    codes = {}
    pattern = os.path.join(harmonization_dir, "taxonomy_cache", "*.json")
    for path in sorted(glob.glob(pattern)):
        code = os.path.splitext(os.path.basename(path))[0]
        with open(path) as fh:
            taxon = json.load(fh).get("taxon", "")
        codes[normalize_binomial(taxon)] = code
    if not codes:
        raise ValueError(f"no taxonomy_cache/*.json under {harmonization_dir} -- is this a stage 02b results dir?")
    return codes


def harmonized_species_codes(harmonization_dir):
    """Every species code that contributed a motif (harmonized_summary.tsv)."""
    path = os.path.join(harmonization_dir, "harmonized_summary.tsv")
    codes = set()
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = header.index("species_code")
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) > idx and fields[idx]:
                codes.add(fields[idx])
    return codes


def resolve_species_codes(manifest, harmonization_dir, allow_unharmonized):
    """{species_id: species_code}, plus a list of warning strings."""
    warnings = []
    if not harmonization_dir:
        return {row["species_id"]: row["species_id"] for row in manifest}, warnings
    cache = load_cache_codes(harmonization_dir)
    contributed = harmonized_species_codes(harmonization_dir)
    resolved = {}
    for row in manifest:
        code = cache.get(normalize_binomial(row["species_name"]))
        if code is None:
            available = sorted(cache)
            msg = (f"species {row['species_id']} ('{row['species_name']}') not found in "
                   f"{harmonization_dir}/taxonomy_cache/*.json (available: {available})")
            if not allow_unharmonized:
                raise ValueError(msg + " -- fix manifest species_name, or set "
                                 "satellite.allow_unharmonized: true to fall back to species_id")
            warnings.append(msg + f"; falling back to species_id '{row['species_id']}'")
            code = row["species_id"]
        elif code not in contributed:
            warnings.append(f"species {row['species_id']} (code {code}) contributed no motifs to "
                            f"harmonized_summary.tsv -- it is masked by the harmonized library "
                            f"without having been part of its discovery")
        resolved[row["species_id"]] = code
    codes = list(resolved.values())
    dupes = sorted({c for c in codes if codes.count(c) > 1})
    if dupes:
        raise ValueError(f"resolved species codes are not unique: {dupes}")
    return resolved, warnings
