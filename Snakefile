# =============================================================================
# Snakefile — repeat_compare: RepeatModeler2 + RepeatMasker cross-sample
# repeat comparison
#
# 1. Sanitizes each manifest sample's genome FASTA, computes assembly QC
#    covariates (contig N50, N content, etc — the denominators used
#    everywhere downstream) and a genome fingerprint (identity guard).
# 2. Runs RepeatModeler2's RECON/RepeatScout rounds (dfam/tetools
#    Singularity image) de novo on each sample independently — WITHOUT
#    -LTRStruct.
# 3. LTR structural discovery outside RepeatModeler: LTR_HARVEST_parallel
#    (+ LTR_FINDER_parallel) on the prepped genome, split into
#    bp-balanced groups of whole scaffolds and fixed-size windows with
#    timeouts, then RepeatModeler's own LTRPipeline downstream steps
#    (LTR_retriever -> MAFFT -> NINJA -> Refiner) on the combined
#    candidates, merged back with the round families exactly the way
#    RepeatModeler does and classified with RepeatClassifier. Only the
#    single whole-genome ltrharvest call (which stalled for days on a
#    highly repetitive hagfish assembly) is replaced.
# 4. Prefixes each sample's family names with its sample_id, optionally
#    exports a Dfam supplement, clusters all samples' families together
#    with cd-hit-est into one non-redundant "shared" library (+ Dfam), and
#    assembles a per-sample "own" library.
# 5. Masks every sample's genome with BOTH libraries (two "arms"), computes
#    divergence landscapes, and summarizes non-overlapping repeat bp per
#    class and per Class/Family against both total and non-N length.
# 6. Combines everything into long-format comparison tables and plots.
#
# A satellite screen / satellite-library arm existed briefly and was removed
# pending fixes to the satellite caller (compare_assemblies_satellites); it
# is preserved at commit 502accf (tag satellite-arm-v1) -- see README.
#
# Sample undergoing programmed germline-to-soma genome rearrangement (e.g.
# hagfish) can have very different repeat content between tissues — this
# pipeline does not act on the manifest's `tissue` column, it only carries
# it into every summary table and warns if samples disagree or are
# `unknown` (see combine_summaries.py).
#
# Usage: snakemake -s Snakefile --configfile config.yaml --cores <N>
# See README.md for full documentation and cluster submission instructions
# (note the --configfile-before-targets caveat documented there).
#
# Author:  KM
# Created: 2026-09
# =============================================================================

import os
import re
import shlex
import sys

configfile: "config.yaml"


# -----------------------------------------------------------------------------
# Manifest parsing (stdlib only, no pandas) — §5.1
# -----------------------------------------------------------------------------
# Header-driven manifest v2: workflow/scripts/manifest.py (shared with
# combine_summaries.py). One row per assembly (sample); see README "Manifest".
sys.path.insert(0, os.path.join(workflow.basedir, "workflow", "scripts"))
from manifest import parse_manifest  # noqa: E402


MANIFEST = parse_manifest(config["manifest"])
SAMPLE_IDS = [row["sample_id"] for row in MANIFEST]
FASTA_BY_SAMPLE = {row["sample_id"]: row["fasta"] for row in MANIFEST}
TISSUE_BY_SAMPLE = {row["sample_id"]: row["tissue"] for row in MANIFEST}
TAXON_BY_SAMPLE = {row["sample_id"]: row["taxon"] for row in MANIFEST}
# Haplotype copies per locus (2 for assembly_type dual_hap): absolute copy
# thresholds are multiplied by it so a two-haplotype assembly isn't judged by
# half the bar.
COPIES_BY_SAMPLE = {row["sample_id"]: row["copies"] for row in MANIFEST}
_dual = [s for s, c in COPIES_BY_SAMPLE.items() if c > 1]
if _dual:
    print(f"[repeat_compare] NOTE: {', '.join(_dual)} are dual_hap assemblies: absolute copy-number "
          f"thresholds (family_tandem.min_copies / major_min_copies / major_min_bp, classify.host_max_copies) are "
          f"doubled for them; percentages are unaffected.")
_sexes = {row["sex"] for row in MANIFEST if row["sex"] != "unknown"}
if len(_sexes) > 1 or any(row["sex"] == "unknown" for row in MANIFEST):
    print(f"[repeat_compare] WARNING: manifest sex values are "
          f"{sorted({row['sex'] for row in MANIFEST})}: W/Y-linked repeats (young ERVs, satellites) "
          f"can't be compared between samples of different or unknown sex.")

tissue_values = set(TISSUE_BY_SAMPLE.values())
if "unknown" in tissue_values or len(tissue_values) > 1:
    print(
        f"[repeat_compare] WARNING: manifest tissue values are {sorted(tissue_values)} — "
        f"a germline assembly and a somatic assembly are different genomes; this is not "
        f"treated as an error, but check before comparing across samples.",
    )

# repeatmasker.sensitive is a single config value used by every arm/sample by
# construction; validate its type here so a future refactor can't silently
# make it per-arm/per-sample without this assertion catching it (§6.8).
if not isinstance(config["repeatmasker"]["sensitive"], bool):
    raise ValueError("config['repeatmasker']['sensitive'] must be true or false")

# Extra RepeatModeler rounds: repeatmodeler.extra_rounds (rule repeatmodeler,
# rm_extend.sh), not -numAddlRounds in extra_args (its meaning changes when a
# run is resumed).
RM_EXTRA_ROUNDS = int(config["repeatmodeler"].get("extra_rounds", 0) or 0)
if "numAddlRounds" in str(config["repeatmodeler"].get("extra_args", "")):
    raise ValueError("repeatmodeler.extra_args: use repeatmodeler.extra_rounds instead of -numAddlRounds "
                     "(extra_rounds: 1 = one more 270 Mb round; an existing run is extended, not redone)")

# -LTRStruct is gone on purpose: LTR discovery runs as the side pipeline
# below (§4 of the restructure plan). Refuse a stale config rather than
# silently ignoring it.
if config["repeatmodeler"].get("ltrstruct"):
    raise ValueError(
        "config repeatmodeler.ltrstruct is no longer supported: RepeatModeler now runs "
        "rounds-only and LTR discovery runs as the ltr_* rules. Remove the key."
    )

OUTDIR = config["outdir"]
ARMS = ["shared", "own"]
INCLUDE_DFAM = bool(config["library"]["include_dfam"])
DFAM_TAXON = config["library"]["dfam_taxon"]
DFAM_EXPORT_FASTA = f"{OUTDIR}/library/dfam_{DFAM_TAXON}.fa"
TETOOLS = config["repeatmodeler"]["container"]

# Mask-only escape hatch: an existing shared library FASTA to mask with
# directly, bypassing discovery/LTR/clustering entirely (rm_library below
# returns it for the shared arm, and no rule produces it, so nothing
# upstream gets scheduled). Meant for the mask_shared_only target, e.g.
#   ./runsnake 40 --configfile config.yaml \
#       --config mask_shared_library=results/library/shared_library.fa \
#       -- mask_shared_only
MASK_SHARED_LIBRARY = config.get("mask_shared_library") or ""
if MASK_SHARED_LIBRARY:
    if not os.path.exists(MASK_SHARED_LIBRARY):
        raise ValueError(f"mask_shared_library does not exist: {MASK_SHARED_LIBRARY}")
    MASK_SHARED_LIBRARY = os.path.abspath(MASK_SHARED_LIBRARY)
    print(
        f"[repeat_compare] WARNING: mask_shared_library is set -- the shared arm masks with "
        f"{MASK_SHARED_LIBRARY} instead of {config['outdir']}/library/shared_library.fa. "
        f"Use it with the mask_shared_only target; unset it for a normal run."
    )

# Prepended to every rule that runs in workflow/envs/repeatmasker.yaml.
# runsnake submits with -V, so a job inherits the PATH of whatever shell
# launched Snakemake -- from a `(base)` shell, miniforge's base python3 can
# come before the rule env's, and famdb.py (`#!/usr/bin/env python3`) then
# dies with "No module named 'h5py'" even though the env has it. Putting the
# rule's own env first also covers everything RepeatMasker calls internally.
# Interpolated via {ENV_PATH_GUARD}; Snakemake doesn't re-format the value,
# so its braces are plain shell.
ENV_PATH_GUARD = """if [ -z "${CONDA_PREFIX:-}" ] || [ ! -x "$CONDA_PREFIX/bin/python" ]; then
    echo "[ERROR] rule conda env not active (CONDA_PREFIX='${CONDA_PREFIX:-}') -- run with --use-conda" >&2
    exit 1
fi
export PATH="$CONDA_PREFIX/bin:$PATH"
echo "[env] CONDA_PREFIX=$CONDA_PREFIX python3=$(command -v python3) famdb.py=$(command -v famdb.py || echo none) RepeatMasker=$(command -v RepeatMasker || echo none)" >&2"""
SCRIPTS = "workflow/scripts"
VENDOR = "workflow/vendor"

# Unknown-family reclassification (README "Reclassifying Unknown families").
# Each screen runs only when its input is configured; an unconfigured screen
# counts as "no evidence" in reclassify_unknown.py.
def _as_bool(value):
    # --config "key={a: false}" delivers nested values as strings ("false"),
    # which bool() would treat as True; config.yaml gives real booleans.
    if isinstance(value, str):
        return value.strip().lower() in ("1", "true", "yes", "on")
    return bool(value)


# Config keys renamed with species -> sample (manifest entries are individual
# assemblies, not species): refuse the old names rather than silently ignore them.
for _sec, _old, _new in (("summary", "plot_species_order", "plot_sample_order"),
                         ("library", "species_prefix_sep", "sample_prefix_sep")):
    if _old in (config.get(_sec) or {}):
        raise ValueError(f"config {_sec}.{_old} was renamed to {_sec}.{_new}; please rename it in config.yaml")

# Left-to-right sample order in the plots (summary.plot_sample_order);
# empty = manifest order. Passed to plot_repeat_compare.R as one argument.
_order = (config.get("summary", {}) or {}).get("plot_sample_order") or []
if isinstance(_order, str):
    _order = [x for x in re.split(r"[,\s\[\]]+", _order) if x]
PLOT_SAMPLE_ORDER = ",".join(list(_order) + [s for s in SAMPLE_IDS if s not in _order])
# Optional groups TSV (family_groups.py members): families that are pieces
# of one element, summed into {summary}/element_groups.tsv.
ELEMENT_GROUPS = config["summary"].get("element_groups", "") or ""
# Optional curated-families table (same format, class_family filled):
# its labels win in summarize, and its elements are grouped too.
CURATED_FAMILIES = (config.get("classify", {}) or {}).get("curated_families", "") or ""
GROUP_TABLES = [p for p in (ELEMENT_GROUPS, CURATED_FAMILIES) if p]

CLASSIFY = config.get("classify", {}) or {}
CLASSIFY_ON = _as_bool(CLASSIFY.get("enabled", False))
ANNOT_PROTEINS = [os.path.abspath(p) for p in (CLASSIFY.get("annotation_proteins") or [])]
SWISSPROT = os.path.abspath(CLASSIFY["swissprot_fasta"]) if CLASSIFY.get("swissprot_fasta") else ""
HOST_SCREEN_ON = CLASSIFY_ON and bool(ANNOT_PROTEINS or SWISSPROT)
RFAM_CM = os.path.abspath(CLASSIFY["rfam_cm"]) if CLASSIFY.get("rfam_cm") else ""
RFAM_ON = CLASSIFY_ON and bool(RFAM_CM)
for _p in ANNOT_PROTEINS + ([SWISSPROT] if SWISSPROT else []) + ([RFAM_CM] if RFAM_CM else []):
    if CLASSIFY_ON and not os.path.exists(_p):
        raise ValueError(f"classify input does not exist: {_p}")


# Satellite / rDNA / mito cross-check (README "Satellite, rDNA and mito
# cross-check"). Optional per-sample inputs come from a header-driven TSV:
# sample_id, tidecluster_dir, tidecluster_prefix, ribotin_fa, mitohifi_fa.
# Every cell must hold a value; ".", "NA", "na", "no", "false", "none" mean
# "not provided" (case-insensitive). Report-only for now (Phase 1).
SATX = config.get("satellite_crosscheck", {}) or {}
SATX_MISSING = {".", "na", "no", "false", "none"}
SATX_COLS = ("sample_id", "tidecluster_dir", "tidecluster_prefix", "ribotin_fa", "mitohifi_fa")


def parse_external_annotations(path):
    rows = {}
    if not path:
        return rows
    header = None
    with open(path) as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split("\t")
            if header is None:
                header = fields
                missing = [c for c in SATX_COLS if c not in header]
                if missing:
                    raise ValueError(f"{path}:{lineno}: header lacks column(s) {missing}; expected {list(SATX_COLS)}")
                continue
            if len(fields) != len(header):
                raise ValueError(f"{path}:{lineno}: {len(fields)} fields, header has {len(header)}")
            for col, val in zip(header, fields):
                if val == "" or val != val.strip() or " " in val:
                    raise ValueError(f"{path}:{lineno}: column '{col}' is blank or contains spaces "
                                     f"(use '.' or NA for 'not provided')")
            row = {c: (None if v.lower() in SATX_MISSING else v) for c, v in zip(header, fields)}
            sid = row["sample_id"]
            if sid not in SAMPLE_IDS:
                raise ValueError(f"{path}:{lineno}: sample_id '{sid}' is not in the manifest")
            if sid in rows:
                raise ValueError(f"{path}:{lineno}: duplicate sample_id '{sid}'")
            if bool(row["tidecluster_dir"]) != bool(row["tidecluster_prefix"]):
                raise ValueError(f"{path}:{lineno}: tidecluster_dir and tidecluster_prefix go together")
            for col in ("tidecluster_dir", "ribotin_fa", "mitohifi_fa"):
                if row[col] and not os.path.exists(row[col]):
                    raise ValueError(f"{path}:{lineno}: {col} does not exist: {row[col]}")
            rows[sid] = row
    return rows


SATX_ANNOT = parse_external_annotations(SATX.get("external_annotations", "") or "")
SATX_ON = bool(SATX_ANNOT)
TC_SAMPLES = [s for s in SAMPLE_IDS if SATX_ANNOT.get(s, {}).get("tidecluster_dir")]
# Samples without TideCluster borrow evidence from a same-taxon sample that has
# it (manifest `taxon`); satellite_crosscheck.infer_from overrides.
SATX_INFER = {s: next(d for d in TC_SAMPLES if TAXON_BY_SAMPLE[d] == TAXON_BY_SAMPLE[s])
              for s in SAMPLE_IDS
              if SATX_ON and s not in TC_SAMPLES and any(TAXON_BY_SAMPLE[d] == TAXON_BY_SAMPLE[s] for d in TC_SAMPLES)}
SATX_INFER.update({str(k): str(v) for k, v in (SATX.get("infer_from") or {}).items()})
for _t, _d in SATX_INFER.items():
    if _t not in SAMPLE_IDS or _d not in TC_SAMPLES:
        raise ValueError(f"satellite_crosscheck.infer_from {_t}: {_d} -- target must be in the manifest "
                         f"and donor must have TideCluster inputs ({TC_SAMPLES})")
if SATX_ON:
    if not TC_SAMPLES:
        raise ValueError("satellite_crosscheck.external_annotations has no sample with TideCluster inputs")
    for _s in SAMPLE_IDS:
        _r = SATX_ANNOT.get(_s, {})
        _mode = ("tidecluster" if _s in TC_SAMPLES
                 else f"inferred from {SATX_INFER[_s]}" if _s in SATX_INFER else "library-level only")
        print(f"[satellite_crosscheck] {_s}: {_mode}; ribotin={'yes' if _r.get('ribotin_fa') else 'no'}; "
              f"mitohifi={'yes' if _r.get('mitohifi_fa') else 'no'}")


def _tc_file(sample, suffix):
    r = SATX_ANNOT[sample]
    return os.path.join(r["tidecluster_dir"], f"{r['tidecluster_prefix']}_{suffix}")


# satellite_crosscheck.apply: summarize also applies the cross-check's
# proposals at apply_confidence (default high), after the user's
# classify.curated_families table (which wins on conflict).
SATX_APPLY = SATX_ON and _as_bool(SATX.get("apply", False))
SATX_APPLIED = f"{OUTDIR}/summary/satellite_applied.tsv"
CURATED_TABLES = ([CURATED_FAMILIES] if CURATED_FAMILIES else []) + ([SATX_APPLIED] if SATX_APPLY else [])

# Optional harmonized satellite motif library (satellite pipeline output):
# names satellites and adds sequence evidence next to the TideCluster arrays.
SATX_MOTIFS = SATX.get("harmonized_library", "") or ""
if SATX_MOTIFS and not os.path.exists(SATX_MOTIFS):
    raise ValueError(f"satellite_crosscheck.harmonized_library not found: {SATX_MOTIFS}")
SATX_REFS = [f"{lab}={SATX_ANNOT[s][col]}" for s in SAMPLE_IDS if s in SATX_ANNOT
             for lab, col in (("rdna", "ribotin_fa"), ("mito", "mitohifi_fa")) if SATX_ANNOT[s][col]]


LTR_CFG = config["ltr_discovery"]
LTR_GROUPS = [f"g{i}" for i in range(int(LTR_CFG["n_groups"]))]
if config["resources"]["ltr_harvest_group"]["threads"] < 2 or config["resources"]["ltr_finder_group"]["threads"] < 2:
    # The vendored *_parallel scripts' 1-thread branch ignores -size/-time
    # (whole group, no timeout) -- exactly the stall this side pipeline avoids.
    raise ValueError("resources.ltr_harvest_group/ltr_finder_group threads must be >= 2")
USE_LTR_FINDER = bool(LTR_CFG["use_ltr_finder"])
LTR_TOOLS = ["harvest", "finder"] if USE_LTR_FINDER else ["harvest"]

wildcard_constraints:
    sample="|".join(re.escape(s) for s in SAMPLE_IDS),
    arm="shared|own",
    group=r"g\d+",
    tool="harvest|finder",

# Chunk count for the repeatmasker scatter/gather split (split_genome /
# repeatmasker_chunk / gather_repeatmasker below) -- reused pattern from
# compare_assemblies_satellites' stage 03 (its -lib-mode RepeatMasker
# scatter/gather), adapted here since our own repeatmasker rule originally
# ran each sample's whole genome as one unchunked job.
scattergather:
    genome_chunks=config["repeatmasker"]["scatter_count"],


# -----------------------------------------------------------------------------
# rule all / library_only — §6.11
# -----------------------------------------------------------------------------
rule all:
    input:
        expand(
            f"{OUTDIR}/{{arm}}/{{sample}}/repeatmasker/{{sample}}.fa.out",
            arm=ARMS,
            sample=SAMPLE_IDS,
        ),
        f"{OUTDIR}/summary/class_composition.tsv",
        f"{OUTDIR}/summary/family_composition.tsv",
        f"{OUTDIR}/summary/divergence_landscape.tsv",
        f"{OUTDIR}/summary/family_tandem.tsv",
        f"{OUTDIR}/summary/class_tandem.tsv",
        f"{OUTDIR}/summary/library_source.tsv",
        [f"{OUTDIR}/summary/element_groups.tsv"] if GROUP_TABLES else [],
        f"{OUTDIR}/summary/curation_candidates.tsv",
        [f"{OUTDIR}/classify/unknown_reclassification.tsv"] if CLASSIFY_ON else [],
        [f"{OUTDIR}/summary/class_verification.tsv", f"{OUTDIR}/summary/class_disagreements.tsv"] if CLASSIFY_ON else [],
        expand(
            f"{OUTDIR}/{{arm}}/{{sample}}/divergence/{{sample}}.landscape.html",
            arm=ARMS,
            sample=SAMPLE_IDS,
        ),
        f"{OUTDIR}/summary/arm_concordance.tsv",
        f"{OUTDIR}/summary/assembly_covariates.tsv",
        f"{OUTDIR}/summary/discovery_round_saturation.tsv",
        f"{OUTDIR}/summary/discovery_round_novelty.tsv",
        f"{OUTDIR}/summary/ltr_discovery.tsv",
        f"{OUTDIR}/summary/ltr_skipped_composition.tsv",
        [f"{OUTDIR}/summary/satellite_family_calls.tsv"] if SATX_ON else [],
        f"{OUTDIR}/summary/provenance.txt",
        f"{OUTDIR}/library/library_membership.tsv",
        f"{OUTDIR}/summary/discovery_summary.tsv",
        [f"{OUTDIR}/summary/dfam_overlap.tsv"] if INCLUDE_DFAM else [],
        f"{OUTDIR}/plots/class_composition_shared.png",
        f"{OUTDIR}/plots/divergence_landscape.png",
        f"{OUTDIR}/plots/arm_concordance.png",


rule library_only:
    # Stops after shared_library.fa so it can be inspected / manually
    # curated (§8) before the expensive masking step.
    input:
        f"{OUTDIR}/library/shared_library.fa",
        f"{OUTDIR}/library/library_membership.tsv",
        f"{OUTDIR}/summary/discovery_summary.tsv",
        [f"{OUTDIR}/summary/dfam_overlap.tsv"] if INCLUDE_DFAM else [],


rule mask_shared_only:
    # RepeatMasker with the shared library on every manifest sample and
    # nothing else (no own arm, no summaries). Pair with mask_shared_library
    # (see top of file) to use an already-built library; without it this
    # still pulls in the full discovery/library chain.
    input:
        expand(
            f"{OUTDIR}/shared/{{sample}}/repeatmasker/{{sample}}.fa.out",
            sample=SAMPLE_IDS,
        ),


def _require_mask_shared_library(target):
    # Input function, so it only fires when the target is actually requested.
    if not MASK_SHARED_LIBRARY:
        raise ValueError(
            f"{target} needs --config mask_shared_library=<the library you masked with> "
            f"(the same value as for mask_shared_only). Without it the shared arm masks "
            f"with {OUTDIR}/library/shared_library.fa, which schedules discovery and "
            f"re-runs RepeatMasker."
        )
    return []


rule report_shared_only:
    # Reporting that needs only the shared arm: divergence + summarize per
    # sample, then combined tables and plots under summary_shared_only/ and
    # plots_shared_only/ (kept apart from a full run's summary/ and plots/).
    # No arm concordance, round saturation, LTR or discovery tables -- those
    # need the own arm or the discovery chain. Requires the same
    # mask_shared_library used for mask_shared_only: without it the shared
    # arm points back at {outdir}/library/shared_library.fa, which pulls in
    # the whole discovery chain AND re-masks every sample.
    input:
        lambda wc: _require_mask_shared_library("report_shared_only"),
        f"{OUTDIR}/summary_shared_only/class_composition.tsv",
        f"{OUTDIR}/summary_shared_only/family_composition.tsv",
        f"{OUTDIR}/summary_shared_only/divergence_landscape.tsv",
        f"{OUTDIR}/summary_shared_only/assembly_covariates.tsv",
        f"{OUTDIR}/summary_shared_only/family_tandem.tsv",
        f"{OUTDIR}/summary_shared_only/class_tandem.tsv",
        f"{OUTDIR}/summary_shared_only/library_source.tsv",
        [f"{OUTDIR}/summary_shared_only/element_groups.tsv"] if GROUP_TABLES else [],
        f"{OUTDIR}/summary_shared_only/curation_candidates.tsv",
        [f"{OUTDIR}/classify/unknown_reclassification.tsv"] if CLASSIFY_ON else [],
        [f"{OUTDIR}/summary_shared_only/class_verification.tsv",
         f"{OUTDIR}/summary_shared_only/class_disagreements.tsv"] if CLASSIFY_ON else [],
        f"{OUTDIR}/plots_shared_only/class_composition_shared.png",
        f"{OUTDIR}/plots_shared_only/divergence_landscape.png",
        expand(
            f"{OUTDIR}/shared/{{sample}}/divergence/{{sample}}.landscape.html",
            sample=SAMPLE_IDS,
        ),


# -----------------------------------------------------------------------------
# 6.1 prep_genome
# -----------------------------------------------------------------------------
rule prep_genome:
    input:
        fasta=lambda wc: FASTA_BY_SAMPLE[wc.sample],
    output:
        fa=f"{OUTDIR}/{{sample}}/genome/{{sample}}.fa",
        name_map=f"{OUTDIR}/{{sample}}/genome/{{sample}}.name_map.tsv",
    threads: config["resources"]["prep_genome"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["prep_genome"]["mem"] * attempt,
        hrs=config["resources"]["prep_genome"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/{{sample}}/prep_genome.log",
    params:
        min_contig_len=config["genome_prep"]["min_contig_len"],
        max_header_len=config["genome_prep"]["max_header_len"],
        test_subsample_bp=config["genome_prep"]["test_subsample_bp"],
    shell:
        "python3 workflow/scripts/prep_genome.py "
        "--fasta {input.fasta} "
        "--min-contig-len {params.min_contig_len} "
        "--max-header-len {params.max_header_len} "
        "--test-subsample-bp {params.test_subsample_bp} "
        "--out-fasta {output.fa} "
        "--out-name-map {output.name_map} "
        "> {log} 2>&1"


# -----------------------------------------------------------------------------
# genome_fingerprint — identity guard (restructure plan: replaces the
# addendum's cross-rule md5 check, which could never differ inside one DAG).
# repeatmodeler compares it against the fingerprint stored when an existing
# RM_* directory was started, before ever -recoverDir-ing it.
# -----------------------------------------------------------------------------
rule genome_fingerprint:
    input:
        fa=f"{OUTDIR}/{{sample}}/genome/{{sample}}.fa",
    output:
        f"{OUTDIR}/{{sample}}/genome/{{sample}}.fingerprint.tsv",
    threads: config["resources"]["genome_fingerprint"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["genome_fingerprint"]["mem"] * attempt,
        hrs=config["resources"]["genome_fingerprint"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/{{sample}}/genome_fingerprint.log",
    shell:
        "python3 {SCRIPTS}/fingerprint.py write --fasta {input.fa} --out {output} > {log} 2>&1"


# -----------------------------------------------------------------------------
# 6.2 assembly_stats
# -----------------------------------------------------------------------------
rule assembly_stats:
    input:
        fa=f"{OUTDIR}/{{sample}}/genome/{{sample}}.fa",
    output:
        f"{OUTDIR}/{{sample}}/genome/{{sample}}.assembly_stats.tsv",
    threads: config["resources"]["assembly_stats"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["assembly_stats"]["mem"] * attempt,
        hrs=config["resources"]["assembly_stats"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/{{sample}}/assembly_stats.log",
    shell:
        "python3 workflow/scripts/assembly_stats.py "
        "--fasta {input.fa} --sample-id {wildcards.sample} --out {output} "
        "> {log} 2>&1"


# -----------------------------------------------------------------------------
# FamDB setup — reused verbatim from compare_assemblies_satellites
# (common/rules/repeatmasker_trf.smk::setup_repeatmasker_famdb), extended to
# also capture RepeatMasker's version and the FamDB release info for
# provenance.txt.
# -----------------------------------------------------------------------------
rule setup_famdb:
    output:
        verified=f"{OUTDIR}/library/famdb_verified.txt",
        rm_version=f"{OUTDIR}/library/repeatmasker_version.txt",
        famdb_release=f"{OUTDIR}/library/famdb_release_info.txt",
    threads: config["resources"]["setup_famdb"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["setup_famdb"]["mem"] * attempt,
        hrs=config["resources"]["setup_famdb"]["hrs"],
        shell_exec="bash",
    conda:
        "workflow/envs/repeatmasker.yaml"
    log:
        f"{OUTDIR}/logs/library/setup_famdb.log",
    params:
        local_dir=config["library"]["famdb_local_dir"].strip(),
        fallback_urls=" ".join(config["library"]["famdb_fallback_urls"]),
        famdb_dir="resources/famdb",
    shell:
        """
        exec > {log} 2>&1
        {ENV_PATH_GUARD}

        RepeatMasker -v > {output.rm_version} 2>&1 || echo "RepeatMasker -v failed" > {output.rm_version}

        write_famdb_conf() {{
            local data_dir="$1"
            while IFS= read -r fpy; do
                fpy_dir="$(dirname "$(readlink -f "$fpy")")"
                cat > "${{fpy_dir}}/famdb.conf" <<CONF
[famdb]
FAMDB_DATA_DIR = $data_dir
CONF
                echo "[INFO] Wrote famdb.conf at ${{fpy_dir}}/famdb.conf -> $data_dir"
            done < <(find "$CONDA_PREFIX" -iname 'famdb.py' 2>/dev/null | sort -u)
        }}

        if famdb.py info > {output.famdb_release} 2>&1; then
            echo "bundled" > {output.verified}
            exit 0
        fi
        echo "[INFO] Bundled 'famdb.py info' failed (expected on a fresh install):"
        cat {output.famdb_release}

        local_dir="{params.local_dir}"
        if [ -n "$local_dir" ]; then
            if famdb.py -i "$local_dir" info > {output.famdb_release} 2>&1; then
                write_famdb_conf "$local_dir"
                echo "local_dir" > {output.verified}
                exit 0
            fi
            echo "[WARN] 'famdb.py -i \"$local_dir\" info' failed — real reason from famdb.py:"
            cat {output.famdb_release}
        fi

        # NOTE: {output.famdb_release} is intentionally left holding whichever
        # attempt above actually ran (local_dir's failure, if one was
        # configured) — it's a diagnostic artifact, not just a release-info
        # cache, so a later failure branch must never blank it.
        fallback_urls="{params.fallback_urls}"
        if [ -z "$fallback_urls" ]; then
            echo "[ERROR] Neither famdb_local_dir nor famdb_fallback_urls is usable." >&2
            echo "[ERROR] See the diagnostic output above (also saved in {output.famdb_release}) for the real reason the local_dir check failed." >&2
            echo "[ERROR] See README.md FamDB section." >&2
            exit 1
        fi

        mkdir -p {params.famdb_dir}
        for url in $fallback_urls; do
            curl -fsSL -o "{params.famdb_dir}/$(basename "$url")" "$url"
        done

        download_dir_abs="$(pwd)/{params.famdb_dir}"
        if famdb.py -i "$download_dir_abs" info > {output.famdb_release} 2>&1; then
            write_famdb_conf "$download_dir_abs"
            echo "fallback_download" > {output.verified}
        else
            echo "[ERROR] famdb.py -i $download_dir_abs info still failed after download." >&2
            exit 1
        fi
        """


# -----------------------------------------------------------------------------
# 6.3 build_db — run inside the SAME Singularity image as repeatmodeler
# (not the conda repeatmasker env), so the BuildDatabase-written database
# files can't drift to a different RepeatModeler/RepeatMasker suite version
# than the one that will actually consume them (spec's §6.3 left this
# rule's environment unspecified; this is the fix applied per the plan).
# -----------------------------------------------------------------------------
rule build_db:
    input:
        fa=f"{OUTDIR}/{{sample}}/genome/{{sample}}.fa",
    output:
        touch(f"{OUTDIR}/{{sample}}/repeatmodeler/{{sample}}.build_db.done"),
    threads: config["resources"]["build_db"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["build_db"]["mem"] * attempt,
        hrs=config["resources"]["build_db"]["hrs"],
        shell_exec="bash",
    singularity:
        TETOOLS
    params:
        workdir=f"{OUTDIR}/{{sample}}/repeatmodeler",
    log:
        f"{OUTDIR}/logs/{{sample}}/build_db.log",
    shell:
        "mkdir -p {params.workdir} && "
        # This container's RepeatModeler 2.0.9 BuildDatabase rejects -engine
        # outright ("Unknown option: engine") -- WU-BLAST support was
        # dropped upstream and NCBI/RMBlast is the only engine now, so the
        # flag is no longer accepted at all, not just unnecessary.
        "BuildDatabase -name {params.workdir}/{wildcards.sample} {input.fa} "
        "> {log} 2>&1"


# -----------------------------------------------------------------------------
# 6.4 repeatmodeler — RECON/RepeatScout rounds only (no -LTRStruct).
#
# Outputs are the rounds' cumulative, UNCLASSIFIED RM_*/consensi.fa and
# families.stk: classification happens once, on the merged rounds + LTR set
# (classify_families), exactly as RepeatModeler itself orders it.
#
# Restart safety:
#  - rm_run.fingerprint.tsv is written next to the RM_* directory right
#    before a fresh run starts. Any later attempt refuses to touch an
#    existing RM_* directory unless the current genome's fingerprint
#    matches it (a -recoverDir on a different genome -- pre-scaffold vs
#    scaffolded, or a changed test_subsample_bp -- would silently mix two
#    genomes).
#  - repeatmodeler.extra_rounds sets the target, 5 + extra_rounds rounds.
#    An interrupted run is recovered in place up to it. A finished run with
#    fewer rounds is copied to RM_*.ext and the copy is extended (the
#    original is never modified; rounds 1-5 must stay identical). See
#    workflow/scripts/rm_extend.sh for RepeatModeler 2.0.9's recovery rules
#    (a round-(h+1)/ dir must exist; resuming at the 270 Mb cap runs
#    -numAddlRounds + 1 rounds).
# -----------------------------------------------------------------------------
rule repeatmodeler:
    input:
        db_done=f"{OUTDIR}/{{sample}}/repeatmodeler/{{sample}}.build_db.done",
        fingerprint=f"{OUTDIR}/{{sample}}/genome/{{sample}}.fingerprint.tsv",
    output:
        consensi=f"{OUTDIR}/{{sample}}/repeatmodeler/{{sample}}.rounds.consensi.fa",
        stk=f"{OUTDIR}/{{sample}}/repeatmodeler/{{sample}}.rounds.families.stk",
        provenance=f"{OUTDIR}/{{sample}}/repeatmodeler/{{sample}}.repeatmodeler_provenance.txt",
    threads: config["resources"]["repeatmodeler"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["repeatmodeler"]["mem"] * attempt,
        hrs=config["resources"]["repeatmodeler"]["hrs"],
        shell_exec="bash",
    singularity:
        TETOOLS
    params:
        workdir=f"{OUTDIR}/{{sample}}/repeatmodeler",
        extra_args=config["repeatmodeler"]["extra_args"],
        extra_rounds=RM_EXTRA_ROUNDS,
        extend_sh=os.path.abspath(f"{SCRIPTS}/rm_extend.sh"),
        fp_abs=lambda wc, input: os.path.abspath(input.fingerprint),
        script_abs=os.path.abspath(f"{SCRIPTS}/fingerprint.py"),
        out_consensi=lambda wc, output: os.path.abspath(output.consensi),
        out_stk=lambda wc, output: os.path.abspath(output.stk),
        out_prov=lambda wc, output: os.path.abspath(output.provenance),
    log:
        f"{OUTDIR}/logs/{{sample}}/repeatmodeler.log",
    shell:
        """
        exec > {log} 2>&1
        set -euo pipefail
        cd {params.workdir}
        bash {params.extend_sh} {wildcards.sample} {threads} {params.extra_rounds} "{params.extra_args}" \
            {params.fp_abs} {params.script_abs} {params.out_consensi} {params.out_stk} {params.out_prov}
        """


# -----------------------------------------------------------------------------
# LTR structural discovery, outside RepeatModeler.
#
# RepeatModeler 2.0.9's -LTRStruct = LTRPipeline: whole-genome gt
# suffixerator + ltrharvest (default parameters, ONE single-threaded
# process -- the step that stalled) -> LTR_retriever -> MAFFT -> NINJA ->
# Refiner. Here only the ltrharvest step is replaced: candidates come from
# LTR_HARVEST_parallel (+ LTR_FINDER_parallel) on the prepped genome, run as bp-balanced groups of whole scaffolds (SGE parallelism)
# and 5 Mb windows with per-window timeouts inside each group (thread
# parallelism). Everything downstream is RepeatModeler's own code
# (workflow/vendor/RepeatModeler/LTRPipeline_from_scn).
# -----------------------------------------------------------------------------
rule ltr_group_genome:
    input:
        f"{OUTDIR}/{{sample}}/genome/{{sample}}.fa",
    output:
        groups=temp(expand(f"{OUTDIR}/{{{{sample}}}}/ltr/groups/{{group}}.fa", group=LTR_GROUPS)),
        manifest=f"{OUTDIR}/{{sample}}/ltr/groups/manifest.tsv",
    threads: config["resources"]["ltr_group_genome"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["ltr_group_genome"]["mem"] * attempt,
        hrs=config["resources"]["ltr_group_genome"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/{{sample}}/ltr_group_genome.log",
    shell:
        "python3 {SCRIPTS}/group_genome.py --fasta {input} --outputs {output.groups} "
        "--manifest {output.manifest} > {log} 2>&1"


_EMPTY_SCN_HEADER = "# no sequences in this group"


rule ltr_harvest_group:
    input:
        f"{OUTDIR}/{{sample}}/ltr/groups/{{group}}.fa",
    output:
        scn=f"{OUTDIR}/{{sample}}/ltr/groups/{{group}}.harvest.scn",
        timeouts=f"{OUTDIR}/{{sample}}/ltr/groups/{{group}}.harvest.timeouts.tsv",
    threads: config["resources"]["ltr_harvest_group"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["ltr_harvest_group"]["mem"] * attempt,
        hrs=config["resources"]["ltr_harvest_group"]["hrs"],
        shell_exec="bash",
    # A stall is handled by the per-window timeouts, not by retrying the
    # whole group with more memory.
    retries: 1
    singularity:
        TETOOLS
    log:
        f"{OUTDIR}/logs/{{sample}}/ltr_harvest_group/{{group}}.log",
    params:
        workdir=f"{OUTDIR}/{{sample}}/ltr/work/harvest_{{group}}",
        fa_abs=lambda wc, input: os.path.abspath(input[0]),
        scn_abs=lambda wc, output: os.path.abspath(output.scn),
        timeouts_abs=lambda wc, output: os.path.abspath(output.timeouts),
        tool=os.path.abspath(f"{VENDOR}/LTR_HARVEST_parallel/LTR_HARVEST_parallel"),
        rm_cfg=os.path.abspath(f"{SCRIPTS}/rm_config_path.sh"),
        size=LTR_CFG["window_size"],
        overlap=LTR_CFG["overlap"],
        time=LTR_CFG["window_timeout_s"],
        try1=LTR_CFG["try1"],
        args=LTR_CFG["ltrharvest_args"],
    shell:
        """
        exec > {log} 2>&1
        set -euo pipefail
        : > {params.timeouts_abs}
        if [ ! -s {params.fa_abs} ]; then
            echo "{_EMPTY_SCN_HEADER}" > {params.scn_abs}
            exit 0
        fi
        # gt is not on the tetools PATH; resolve it the way RepeatModeler does.
        GT_DIR=$(bash {params.rm_cfg} GENOMETOOLS_DIR)
        rm -rf {params.workdir} && mkdir -p {params.workdir} && cd {params.workdir}
        perl {params.tool} -seq {params.fa_abs} -size {params.size} -overlap {params.overlap} \
            -time {params.time} -try1 {params.try1} -threads {threads} -gt "$GT_DIR" \
            -harvest_args "{params.args}" -timeout_log {params.timeouts_abs}
        cp {wildcards.group}.fa.harvest.combine.scn {params.scn_abs}
        cd - > /dev/null && rm -rf {params.workdir}
        """


rule ltr_finder_group:
    input:
        f"{OUTDIR}/{{sample}}/ltr/groups/{{group}}.fa",
    output:
        scn=f"{OUTDIR}/{{sample}}/ltr/groups/{{group}}.finder.scn",
        timeouts=f"{OUTDIR}/{{sample}}/ltr/groups/{{group}}.finder.timeouts.tsv",
        version=f"{OUTDIR}/{{sample}}/ltr/groups/{{group}}.finder.version.txt",
    threads: config["resources"]["ltr_finder_group"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["ltr_finder_group"]["mem"] * attempt,
        hrs=config["resources"]["ltr_finder_group"]["hrs"],
        shell_exec="bash",
    retries: 1
    conda:
        "workflow/envs/ltr_finder.yaml"
    log:
        f"{OUTDIR}/logs/{{sample}}/ltr_finder_group/{{group}}.log",
    params:
        workdir=f"{OUTDIR}/{{sample}}/ltr/work/finder_{{group}}",
        fa_abs=lambda wc, input: os.path.abspath(input[0]),
        scn_abs=lambda wc, output: os.path.abspath(output.scn),
        timeouts_abs=lambda wc, output: os.path.abspath(output.timeouts),
        version_abs=lambda wc, output: os.path.abspath(output.version),
        tool=os.path.abspath(f"{VENDOR}/LTR_FINDER_parallel/LTR_FINDER_parallel"),
        size=LTR_CFG["window_size"],
        overlap=LTR_CFG["overlap"],
        time=LTR_CFG["window_timeout_s"],
        try1=LTR_CFG["try1"],
        args=LTR_CFG["ltr_finder_args"],
    shell:
        """
        exec > {log} 2>&1
        set -euo pipefail
        : > {params.timeouts_abs}
        (ltr_finder 2>&1 | grep -i -m1 version || echo "ltr_finder (version line not found)") > {params.version_abs}
        if [ ! -s {params.fa_abs} ]; then
            echo "{_EMPTY_SCN_HEADER}" > {params.scn_abs}
            exit 0
        fi
        perl -Mthreads -e 1 || {{ echo "[ERROR] this perl lacks ithreads; LTR_FINDER_parallel needs them"; exit 1; }}
        FINDER_DIR=$(dirname "$(readlink -f "$(command -v ltr_finder)")")
        rm -rf {params.workdir} && mkdir -p {params.workdir} && cd {params.workdir}
        perl {params.tool} -seq {params.fa_abs} -size {params.size} -overlap {params.overlap} \
            -time {params.time} -try1 {params.try1} -threads {threads} -harvest_out \
            -finder "$FINDER_DIR" -finder_args "{params.args}" -timeout_log {params.timeouts_abs}
        cp {wildcards.group}.fa.finder.combine.scn {params.scn_abs}
        cd - > /dev/null && rm -rf {params.workdir}
        """


rule ltr_gather:
    input:
        genome=f"{OUTDIR}/{{sample}}/genome/{{sample}}.fa",
        harvest=expand(f"{OUTDIR}/{{{{sample}}}}/ltr/groups/{{group}}.harvest.scn", group=LTR_GROUPS),
        harvest_to=expand(f"{OUTDIR}/{{{{sample}}}}/ltr/groups/{{group}}.harvest.timeouts.tsv", group=LTR_GROUPS),
        finder=expand(f"{OUTDIR}/{{{{sample}}}}/ltr/groups/{{group}}.finder.scn", group=LTR_GROUPS) if USE_LTR_FINDER else [],
        finder_to=expand(f"{OUTDIR}/{{{{sample}}}}/ltr/groups/{{group}}.finder.timeouts.tsv", group=LTR_GROUPS) if USE_LTR_FINDER else [],
    output:
        scn=f"{OUTDIR}/{{sample}}/ltr/rawLTR.scn",
        skipped=f"{OUTDIR}/{{sample}}/ltr/skipped_windows.tsv",
        summary=f"{OUTDIR}/{{sample}}/ltr/ltr_discovery_summary.tsv",
    threads: config["resources"]["ltr_gather"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["ltr_gather"]["mem"] * attempt,
        hrs=config["resources"]["ltr_gather"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/{{sample}}/ltr_gather.log",
    params:
        timeout_logs=lambda wc, input: " ".join(
            [f"harvest:{p}" for p in input.harvest_to] + [f"finder:{p}" for p in input.finder_to]
        ),
        finder_arg=lambda wc, input: f"--finder {' '.join(input.finder)}" if input.finder else "",
        size=LTR_CFG["window_size"],
        overlap=LTR_CFG["overlap"],
    shell:
        "python3 {SCRIPTS}/normalize_scn.py --genome {input.genome} --harvest {input.harvest} "
        "{params.finder_arg} --timeout-logs {params.timeout_logs} "
        "--window-size {params.size} --overlap {params.overlap} --sample {wildcards.sample} "
        "--out-scn {output.scn} --out-skipped {output.skipped} --out-summary {output.summary} "
        "> {log} 2>&1"


rule ltr_pipeline:
    input:
        genome=f"{OUTDIR}/{{sample}}/genome/{{sample}}.fa",
        scn=f"{OUTDIR}/{{sample}}/ltr/rawLTR.scn",
    output:
        fa=f"{OUTDIR}/{{sample}}/ltr/{{sample}}.ltrs.fa",
        stk=f"{OUTDIR}/{{sample}}/ltr/{{sample}}.ltrs.stk",
        versions=f"{OUTDIR}/{{sample}}/ltr/tool_versions.txt",
    threads: config["resources"]["ltr_pipeline"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["ltr_pipeline"]["mem"] * attempt,
        hrs=config["resources"]["ltr_pipeline"]["hrs"],
        shell_exec="bash",
    singularity:
        TETOOLS
    log:
        f"{OUTDIR}/logs/{{sample}}/ltr_pipeline.log",
    params:
        workdir=f"{OUTDIR}/{{sample}}/ltr/work/pipeline",
        genome_abs=lambda wc, input: os.path.abspath(input.genome),
        scn_abs=lambda wc, input: os.path.abspath(input.scn),
        fa_abs=lambda wc, output: os.path.abspath(output.fa),
        stk_abs=lambda wc, output: os.path.abspath(output.stk),
        versions_abs=lambda wc, output: os.path.abspath(output.versions),
        tool=os.path.abspath(f"{VENDOR}/RepeatModeler/LTRPipeline_from_scn"),
        rm_cfg=os.path.abspath(f"{SCRIPTS}/rm_config_path.sh"),
    shell:
        """
        exec > {log} 2>&1
        set -euo pipefail
        {{
            echo "genometools: $("$(bash {params.rm_cfg} GENOMETOOLS_DIR)/gt" --version 2>&1 | head -n1)"
            echo "LTR_retriever: $(bash {params.rm_cfg} LTR_RETRIEVER_DIR)"
            grep -m1 -i "version" "$(bash {params.rm_cfg} LTR_RETRIEVER_DIR)/LTR_retriever" || true
            echo "RepeatModeler: $(RepeatModeler -version 2>&1 | head -n1)"
        }} > {params.versions_abs}
        rm -rf {params.workdir} && mkdir -p {params.workdir} && cd {params.workdir}
        # LTRPipeline writes <input>-ltrs.fa next to its input; link the
        # genome in so everything stays inside the work dir.
        ln -s {params.genome_abs} {wildcards.sample}.fa
        export RM_DIR=$(dirname "$(readlink -f "$(command -v RepeatModeler)")")
        perl {params.tool} -inscn {params.scn_abs} -threads {threads} -tmpdir . {wildcards.sample}.fa
        if [ -s {wildcards.sample}.fa-ltrs.fa ]; then
            cp {wildcards.sample}.fa-ltrs.fa {params.fa_abs}
            cp {wildcards.sample}.fa-ltrs.stk {params.stk_abs}
        else
            echo "[WARN] LTRPipeline produced no LTR families (see above); continuing rounds-only"
            : > {params.fa_abs}
            : > {params.stk_abs}
        fi
        cd - > /dev/null && rm -rf {params.workdir}
        """


# -----------------------------------------------------------------------------
# Merge back + classify, as RepeatModeler does after -LTRStruct
# -----------------------------------------------------------------------------
rule merge_families:
    input:
        rounds_fa=f"{OUTDIR}/{{sample}}/repeatmodeler/{{sample}}.rounds.consensi.fa",
        rounds_stk=f"{OUTDIR}/{{sample}}/repeatmodeler/{{sample}}.rounds.families.stk",
        ltr_fa=f"{OUTDIR}/{{sample}}/ltr/{{sample}}.ltrs.fa",
        ltr_stk=f"{OUTDIR}/{{sample}}/ltr/{{sample}}.ltrs.stk",
    output:
        fa=f"{OUTDIR}/{{sample}}/families/{{sample}}.merged.consensi.fa",
        stk=f"{OUTDIR}/{{sample}}/families/{{sample}}.merged.families.stk",
    threads: config["resources"]["merge_families"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["merge_families"]["mem"] * attempt,
        hrs=config["resources"]["merge_families"]["hrs"],
        shell_exec="bash",
    singularity:
        TETOOLS
    log:
        f"{OUTDIR}/logs/{{sample}}/merge_families.log",
    params:
        workdir=f"{OUTDIR}/{{sample}}/families/merge_work",
        rm_cfg=f"{SCRIPTS}/rm_config_path.sh",
    shell:
        """
        exec > {log} 2>&1
        set -euo pipefail
        CDHIT="$(bash {params.rm_cfg} CDHIT_DIR)/cd-hit-est"
        python3 {SCRIPTS}/merge_families.py --rounds-fa {input.rounds_fa} --rounds-stk {input.rounds_stk} \
            --ltr-fa {input.ltr_fa} --ltr-stk {input.ltr_stk} --cdhit "$CDHIT" --threads {threads} \
            --workdir {params.workdir} --out-fa {output.fa} --out-stk {output.stk}
        rm -rf {params.workdir}
        """


rule classify_families:
    input:
        fa=f"{OUTDIR}/{{sample}}/families/{{sample}}.merged.consensi.fa",
        stk=f"{OUTDIR}/{{sample}}/families/{{sample}}.merged.families.stk",
    output:
        fa=f"{OUTDIR}/{{sample}}/families/{{sample}}-families.fa",
        stk=f"{OUTDIR}/{{sample}}/families/{{sample}}-families.stk",
    threads: config["resources"]["classify_families"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["classify_families"]["mem"] * attempt,
        hrs=config["resources"]["classify_families"]["hrs"],
        shell_exec="bash",
    singularity:
        TETOOLS
    log:
        f"{OUTDIR}/logs/{{sample}}/classify_families.log",
    params:
        workdir=f"{OUTDIR}/{{sample}}/families/classify_work",
        fa_abs=lambda wc, output: os.path.abspath(output.fa),
        stk_abs=lambda wc, output: os.path.abspath(output.stk),
    shell:
        """
        exec > {log} 2>&1
        set -euo pipefail
        rm -rf {params.workdir} && mkdir -p {params.workdir}
        cp {input.fa} {params.workdir}/consensi.fa
        cp {input.stk} {params.workdir}/families.stk
        cd {params.workdir}
        RepeatClassifier -consensi consensi.fa -stockholm families.stk -threads {threads}
        cp consensi.fa.classified {params.fa_abs}
        cp families-classified.stk {params.stk_abs}
        cd - > /dev/null && rm -rf {params.workdir}
        """


# -----------------------------------------------------------------------------
# 6.5 prefix_library
# -----------------------------------------------------------------------------
rule prefix_library:
    input:
        fa=f"{OUTDIR}/{{sample}}/families/{{sample}}-families.fa",
    output:
        f"{OUTDIR}/{{sample}}/library/{{sample}}.prefixed.fa",
    threads: config["resources"]["prefix_library"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["prefix_library"]["mem"] * attempt,
        hrs=config["resources"]["prefix_library"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/{{sample}}/prefix_library.log",
    params:
        sep=config["library"]["sample_prefix_sep"],
    shell:
        "python3 workflow/scripts/prefix_library.py "
        "--fasta {input.fa} --sample-id {wildcards.sample} --sep {params.sep} "
        "--out {output} > {log} 2>&1"


# -----------------------------------------------------------------------------
# 6.6 export_dfam (only if library.include_dfam)
# -----------------------------------------------------------------------------
if INCLUDE_DFAM:

    rule export_dfam:
        input:
            verified=f"{OUTDIR}/library/famdb_verified.txt",
        output:
            DFAM_EXPORT_FASTA,
        threads: config["resources"]["export_dfam"]["threads"]
        resources:
            mem=lambda wildcards, attempt: config["resources"]["export_dfam"]["mem"] * attempt,
            hrs=config["resources"]["export_dfam"]["hrs"],
            shell_exec="bash",
        conda:
            "workflow/envs/repeatmasker.yaml"
        log:
            f"{OUTDIR}/logs/library/export_dfam.log",
        params:
            taxon=DFAM_TAXON,
        shell:
            # -f fasta is not a valid --format choice on this installed
            # famdb.py (confirmed via the actual usage error: choices are
            # summary/hmm/hmm_species/fasta_name/fasta_acc/embl*) --
            # fasta_name gives human-readable family-name headers matching
            # the rest of this pipeline's library naming, and
            # --include-class-in-name (a separate flag from --format) is
            # what actually appends the #Class/Family suffix (spec §6.6).
            #
            # -c/--curated (confirmed via `famdb.py families -h`): without
            # it, -a -d on a broad taxon like Vertebrata returns EVERY
            # uncurated per-genome RepeatModeler-derived family Dfam has
            # ever ingested for that clade -- 4.4M entries in practice, not
            # a "supplement" (spec §6.6) by any reading. -c restricts to
            # the curated cross-sample reference set that section actually
            # describes.
            #
            # No --add-reverse-complement: it writes every family twice
            # (forward + "(anti)"), and RepeatMasker already searches both
            # strands, so the copies only doubled the Dfam search cost.
            """
            exec 2> {log}
            {ENV_PATH_GUARD}
            famdb.py families -f fasta_name --include-class-in-name -c -a -d \
                '{params.taxon}' > {output}
            """


# -----------------------------------------------------------------------------
# 6.7 cluster_library + library_membership + shared/own library assembly
# -----------------------------------------------------------------------------
rule cluster_library:
    input:
        prefixed=expand(f"{OUTDIR}/{{sample}}/library/{{sample}}.prefixed.fa", sample=SAMPLE_IDS),
    output:
        all_prefixed=temp(f"{OUTDIR}/library/all_prefixed.fa"),
        nr_fa=f"{OUTDIR}/library/shared_denovo.nr.fa",
        clstr=f"{OUTDIR}/library/shared_denovo.nr.fa.clstr",
        cdhit_version=f"{OUTDIR}/library/cdhit_version.txt",
    threads: config["resources"]["cluster_library"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["cluster_library"]["mem"] * attempt,
        hrs=config["resources"]["cluster_library"]["hrs"],
        shell_exec="bash",
    conda:
        "workflow/envs/cdhit.yaml"
    log:
        f"{OUTDIR}/logs/library/cluster_library.log",
    params:
        identity=config["library"]["cdhit"]["identity"],
        coverage_short=config["library"]["cdhit"]["coverage_short"],
        word_size=config["library"]["cdhit"]["word_size"],
        total_mb=lambda wc, threads, resources: resources.mem * threads * 1024,
    shell:
        "cat {input.prefixed} > {output.all_prefixed} && "
        "cd-hit-est -i {output.all_prefixed} -o {output.nr_fa} "
        "-c {params.identity} -aS {params.coverage_short} -n {params.word_size} "
        "-G 0 -d 0 -r 1 -M {params.total_mb} -T {threads} "
        "> {log} 2>&1 && "
        "(cd-hit-est 2>&1 | head -n1) > {output.cdhit_version} || true"


rule library_membership:
    input:
        clstr=f"{OUTDIR}/library/shared_denovo.nr.fa.clstr",
        dfam=[DFAM_EXPORT_FASTA] if INCLUDE_DFAM else [],
    output:
        f"{OUTDIR}/library/library_membership.tsv",
    threads: config["resources"]["library_membership"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["library_membership"]["mem"] * attempt,
        hrs=config["resources"]["library_membership"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/library/library_membership.log",
    params:
        sep=config["library"]["sample_prefix_sep"],
        codes=" ".join(SAMPLE_IDS),
        dfam_arg=lambda wc, input: f"--dfam {input.dfam}" if input.dfam else "",
    shell:
        "python3 workflow/scripts/library_membership.py "
        "--clstr {input.clstr} --sep {params.sep} --sample-codes {params.codes} "
        "{params.dfam_arg} --out {output} > {log} 2>&1"


rule discovery_summary:
    # Family counts per sample from RepeatModeler rounds + LTR pipeline
    # through merge, classification and cross-sample clustering (shared
    # vs sample-only), plus a per-Class breakdown. Library-stage inputs
    # only, so it's also built by library_only.
    input:
        rounds_fa=expand(f"{OUTDIR}/{{sample}}/repeatmodeler/{{sample}}.rounds.consensi.fa", sample=SAMPLE_IDS),
        ltr_fa=expand(f"{OUTDIR}/{{sample}}/ltr/{{sample}}.ltrs.fa", sample=SAMPLE_IDS),
        merged_fa=expand(f"{OUTDIR}/{{sample}}/families/{{sample}}.merged.consensi.fa", sample=SAMPLE_IDS),
        classified_fa=expand(f"{OUTDIR}/{{sample}}/families/{{sample}}-families.fa", sample=SAMPLE_IDS),
        clstr=f"{OUTDIR}/library/shared_denovo.nr.fa.clstr",
        membership=f"{OUTDIR}/library/library_membership.tsv",
    output:
        summary=f"{OUTDIR}/summary/discovery_summary.tsv",
        by_class=f"{OUTDIR}/summary/discovery_summary_by_class.tsv",
    threads: config["resources"]["discovery_summary"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["discovery_summary"]["mem"] * attempt,
        hrs=config["resources"]["discovery_summary"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/summary/discovery_summary.log",
    params:
        sep=config["library"]["sample_prefix_sep"],
        sample=" ".join(SAMPLE_IDS),
        taxon=" ".join(f"{s}={TAXON_BY_SAMPLE[s]}" for s in SAMPLE_IDS),
    shell:
        "python3 {SCRIPTS}/discovery_summary.py --sample {params.sample} --taxon {params.taxon} "
        "--rounds-fa {input.rounds_fa} --ltr-fa {input.ltr_fa} --merged-fa {input.merged_fa} "
        "--classified-fa {input.classified_fa} --clstr {input.clstr} --membership {input.membership} "
        "--sep {params.sep} --out {output.summary} --by-class-out {output.by_class} > {log} 2>&1"


if INCLUDE_DFAM:

    rule dfam_overlap_sample:
        # Which of this sample's de novo families are already-known Dfam
        # families: cd-hit-est-2d with db1 = the Dfam export, db2 = the
        # sample's prefixed families, same thresholds as cluster_library.
        # -s2 0.8: a Dfam family must be >= 80% of the de novo family's
        # length. -aS alone (coverage of the SHORTER sequence) let any family
        # that merely contains a short Dfam entry (tRNA, MITE, solo LTR)
        # match it; -aL doesn't help in cd-hit-2d. -S2 lifts the default
        # "db2 no longer than db1" so a family with some flank still joins.
        # dfam_overlap.py then requires >= 80% coverage of the de novo family.
        input:
            dfam=DFAM_EXPORT_FASTA,
            families=f"{OUTDIR}/{{sample}}/library/{{sample}}.prefixed.fa",
        output:
            clstr=f"{OUTDIR}/library/dfam_overlap/{{sample}}.dfam2d.clstr",
            unmatched=temp(f"{OUTDIR}/library/dfam_overlap/{{sample}}.dfam2d"),
        threads: config["resources"]["dfam_overlap"]["threads"]
        resources:
            mem=lambda wildcards, attempt: config["resources"]["dfam_overlap"]["mem"] * attempt,
            hrs=config["resources"]["dfam_overlap"]["hrs"],
            shell_exec="bash",
        conda:
            "workflow/envs/cdhit.yaml"
        log:
            f"{OUTDIR}/logs/library/dfam_overlap_{{sample}}.log",
        params:
            identity=config["library"]["cdhit"]["identity"],
            coverage_short=config["library"]["cdhit"]["coverage_short"],
            word_size=config["library"]["cdhit"]["word_size"],
            total_mb=lambda wc, threads, resources: resources.mem * threads * 1024,
        shell:
            "cd-hit-est-2d -i {input.dfam} -i2 {input.families} -o {output.unmatched} "
            "-c {params.identity} -aS {params.coverage_short} -n {params.word_size} "
            "-G 0 -g 1 -r 1 -d 0 -s2 0.8 -S2 999999999 -M {params.total_mb} -T {threads} "
            "> {log} 2>&1"

    rule dfam_overlap:
        input:
            families=expand(f"{OUTDIR}/{{sample}}/library/{{sample}}.prefixed.fa", sample=SAMPLE_IDS),
            clstr=expand(f"{OUTDIR}/library/dfam_overlap/{{sample}}.dfam2d.clstr", sample=SAMPLE_IDS),
        output:
            summary=f"{OUTDIR}/summary/dfam_overlap.tsv",
            matches=f"{OUTDIR}/library/dfam_overlap/dfam_matches.tsv",
        threads: config["resources"]["discovery_summary"]["threads"]
        resources:
            mem=lambda wildcards, attempt: config["resources"]["discovery_summary"]["mem"] * attempt,
            hrs=config["resources"]["discovery_summary"]["hrs"],
            shell_exec="bash",
        log:
            f"{OUTDIR}/logs/summary/dfam_overlap.log",
        params:
            sample=" ".join(SAMPLE_IDS),
            min_coverage=config["library"]["cdhit"]["coverage_short"],
        shell:
            "python3 {SCRIPTS}/dfam_overlap.py --sample {params.sample} "
            "--families {input.families} --clstr {input.clstr} "
            "--min-coverage {params.min_coverage} "
            "--out-summary {output.summary} --out-matches {output.matches} > {log} 2>&1"


def _shared_library_inputs(wildcards):
    inputs = {"nr": f"{OUTDIR}/library/shared_denovo.nr.fa"}
    if INCLUDE_DFAM:
        inputs["dfam"] = DFAM_EXPORT_FASTA
    return inputs


rule assemble_shared_library:
    # Primary comparison-arm library: clustered de novo families (or
    # library.curated_override, which replaces clustering for THIS file,
    # §8) + the optional Dfam export. cluster_library and
    # library_membership.tsv still run unconditionally, since the de novo
    # shared-vocabulary comparison is a result in its own right.
    input:
        unpack(_shared_library_inputs),
    output:
        fa=f"{OUTDIR}/library/shared_library.fa",
        report=f"{OUTDIR}/library/shared_library.sources.tsv",
    threads: config["resources"]["assemble_shared_library"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["assemble_shared_library"]["mem"] * attempt,
        hrs=config["resources"]["assemble_shared_library"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/library/assemble_shared_library.log",
    params:
        base=lambda wc, input: config["library"]["curated_override"] or input.nr,
        dfam_arg=lambda wc, input: f"--dfam {input.dfam}" if INCLUDE_DFAM else "",
    shell:
        "python3 {SCRIPTS}/append_libraries.py --base {params.base} {params.dfam_arg} "
        "--out {output.fa} --report {output.report} > {log} 2>&1"


def _own_library_inputs(wildcards):
    inputs = {"prefixed": f"{OUTDIR}/{wildcards.sample}/library/{wildcards.sample}.prefixed.fa"}
    if INCLUDE_DFAM:
        inputs["dfam"] = DFAM_EXPORT_FASTA
    return inputs


rule own_library:
    # Sanity-check arm library: that sample's de novo families only (+ the
    # same optional Dfam export), assembled in its own small rule so both
    # arms call RepeatMasker identically (§6.7).
    input:
        unpack(_own_library_inputs),
    output:
        fa=f"{OUTDIR}/own/{{sample}}/library/{{sample}}.own_library.fa",
        report=f"{OUTDIR}/own/{{sample}}/library/{{sample}}.own_library.sources.tsv",
    threads: config["resources"]["own_library"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["own_library"]["mem"] * attempt,
        hrs=config["resources"]["own_library"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/{{sample}}/own_library.log",
    params:
        dfam_arg=lambda wc, input: f"--dfam {input.dfam}" if INCLUDE_DFAM else "",
    shell:
        "python3 {SCRIPTS}/append_libraries.py --base {input.prefixed} {params.dfam_arg} "
        "--out {output.fa} --report {output.report} > {log} 2>&1"


# -----------------------------------------------------------------------------
# 6.8 repeatmasker (both arms) -- scatter/gather over genome chunks
#
# Adapted from compare_assemblies_satellites' stage 03 -lib-mode RepeatMasker
# scatter/gather (split_genome_fasta / run_repeatmasker_genome_lib /
# gather_repeatmasker_genome_lib), rather than running each sample's whole
# genome as one unchunked job: real multi-Gb assemblies schedule much better
# on SGE as N independently-restartable chunk jobs. Deliberately NOT copying
# that stage's -nolow flag -- it suppresses RepeatMasker's built-in
# low-complexity/simple-repeat screen, which we need (Simple_repeat and
# Low_complexity are both required canonical output classes here).
# -----------------------------------------------------------------------------
def rm_library(wildcards):
    if wildcards.arm == "shared":
        return MASK_SHARED_LIBRARY or f"{OUTDIR}/library/shared_library.fa"
    return f"{OUTDIR}/own/{wildcards.sample}/library/{wildcards.sample}.own_library.fa"


rule split_genome:
    # Scatters one sample's genome into config["repeatmasker"]["scatter_count"]
    # chunks (workflow/scripts/split_fasta.py, copied verbatim from the
    # sibling repo's common/scripts/split_fasta.py -- snake/boustrophedon by
    # contig count). Per-sample, not per-arm: the genome being split is
    # identical regardless of which library later masks it.
    input:
        fasta=f"{OUTDIR}/{{sample}}/genome/{{sample}}.fa",
    output:
        fasta=temp(scatter.genome_chunks(
            f"{OUTDIR}/{{{{sample}}}}/genome/chunks/{{scatteritem}}/{{scatteritem}}.fa"
        )),
    threads: config["resources"]["split_genome"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["split_genome"]["mem"] * attempt,
        hrs=config["resources"]["split_genome"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/{{sample}}/split_genome.log",
    shell:
        "python3 workflow/scripts/split_fasta.py --infile {input.fasta} "
        "--outputs {output.fasta} > {log} 2>&1"


rule repeatmasker_chunk:
    input:
        fasta=f"{OUTDIR}/{{sample}}/genome/chunks/{{scatteritem}}/{{scatteritem}}.fa",
        lib=rm_library,
        famdb_verified=f"{OUTDIR}/library/famdb_verified.txt",
    output:
        out_file=temp(f"{OUTDIR}/{{arm}}/{{sample}}/repeatmasker/chunks/{{scatteritem}}/{{scatteritem}}.fa.out"),
        tbl_file=temp(f"{OUTDIR}/{{arm}}/{{sample}}/repeatmasker/chunks/{{scatteritem}}/{{scatteritem}}.fa.tbl"),
        align_file=temp(f"{OUTDIR}/{{arm}}/{{sample}}/repeatmasker/chunks/{{scatteritem}}/{{scatteritem}}.fa.align"),
    threads: config["resources"]["repeatmasker"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["repeatmasker"]["mem"] * attempt,
        hrs=config["resources"]["repeatmasker"]["hrs"],
        shell_exec="bash",
    conda:
        "workflow/envs/repeatmasker.yaml"
    params:
        outdir=f"{OUTDIR}/{{arm}}/{{sample}}/repeatmasker/chunks/{{scatteritem}}",
        pa=config["resources"]["repeatmasker"]["threads"] // config["repeatmasker"]["cores_per_pa"],
        sensitive_flag="-s" if config["repeatmasker"]["sensitive"] else "",
        extra_args=config["repeatmasker"]["extra_args"],
        # Resolved here (DAG-build time, original invocation directory) --
        # NOT via a shell-level `readlink -f` inside the rule body, which
        # would run after the `cd {params.outdir}` below and resolve these
        # relative paths against the wrong directory, silently expanding to
        # nothing (readlink -f fails when the leading path components don't
        # exist, and a failed command substitution doesn't trip `set -e`).
        fa_abs=lambda wc, input: os.path.abspath(input.fasta),
        lib_abs=lambda wc, input: os.path.abspath(input.lib),
    log:
        f"{OUTDIR}/logs/{{arm}}/{{sample}}/repeatmasker_chunk/{{scatteritem}}.log",
    shell:
        """
        exec > {log} 2>&1
        {ENV_PATH_GUARD}
        mkdir -p {params.outdir}
        (cd {params.outdir} && trap 'rm -rf RM_*' EXIT && \
            RepeatMasker -pa {params.pa} -lib {params.lib_abs} -xsmall -gff -a \
                {params.sensitive_flag} {params.extra_args} -dir . {params.fa_abs})
        touch {output.out_file} {output.tbl_file} {output.align_file}
        """


rule gather_repeatmasker:
    # .out: workflow/scripts/gather_rm_out.sh, copied verbatim from the
    # sibling repo's common/scripts/gather_rm_out.sh -- keeps one chunk's
    # 3-line header, appends every chunk's data rows, strips the
    # "There were no repetitive sequences detected" zero-hit sentinel line
    # per chunk (not just once), so it merges cleanly regardless of which
    # chunks had zero hits.
    # .tbl: the sibling repo never needed this (stage 03 only merges .out).
    # Chunks are disjoint contig sets, so "total length"/"bases masked" are
    # additive -- sum them across chunks into a minimal synthetic .tbl
    # (summarize_rm.py's parse_tbl_masked_bp() only greps "bases masked"
    # out of it anyway, so a two-line file is sufficient).
    # .align: also new -- straight concatenation (no shared header block
    # like .out has). Verify this against real wiring-test output before
    # trusting it for the full run.
    input:
        out_chunks=gather.genome_chunks(
            f"{OUTDIR}/{{{{arm}}}}/{{{{sample}}}}/repeatmasker/chunks/{{scatteritem}}/{{scatteritem}}.fa.out"
        ),
        tbl_chunks=gather.genome_chunks(
            f"{OUTDIR}/{{{{arm}}}}/{{{{sample}}}}/repeatmasker/chunks/{{scatteritem}}/{{scatteritem}}.fa.tbl"
        ),
        align_chunks=gather.genome_chunks(
            f"{OUTDIR}/{{{{arm}}}}/{{{{sample}}}}/repeatmasker/chunks/{{scatteritem}}/{{scatteritem}}.fa.align"
        ),
    output:
        out_file=f"{OUTDIR}/{{arm}}/{{sample}}/repeatmasker/{{sample}}.fa.out",
        tbl_file=f"{OUTDIR}/{{arm}}/{{sample}}/repeatmasker/{{sample}}.fa.tbl",
        align_file=f"{OUTDIR}/{{arm}}/{{sample}}/repeatmasker/{{sample}}.fa.align",
    threads: config["resources"]["gather_repeatmasker"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["gather_repeatmasker"]["mem"] * attempt,
        hrs=config["resources"]["gather_repeatmasker"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/{{arm}}/{{sample}}/gather_repeatmasker.log",
    shell:
        """
        exec > {log} 2>&1
        bash workflow/scripts/gather_rm_out.sh {output.out_file} {input.out_chunks}
        cat {input.align_chunks} > {output.align_file}
        total_len=$(awk -F'[ :]+' '/total length/{{sum+=$3}} END{{print sum+0}}' {input.tbl_chunks})
        masked_bp=$(awk -F'[ :]+' '/bases masked/{{sum+=$3}} END{{print sum+0}}' {input.tbl_chunks})
        {{
            echo "total length: ${{total_len}} bp"
            echo "bases masked: ${{masked_bp}} bp"
        }} > {output.tbl_file}
        """


# -----------------------------------------------------------------------------
# 6.9 divergence
#
# RepeatMasker's stock .divsum / landscape .html, kept for reference only.
# calcDivergenceFromAlign.pl sums every alignment in .align, including the
# overlapping ones .out resolves (tandem arrays especially), so its bp can
# exceed the genome. summary/divergence_landscape.tsv and the plot come from
# summarize_rm.py's overlap-resolved landscape instead.
# -----------------------------------------------------------------------------
rule divergence:
    input:
        align=f"{OUTDIR}/{{arm}}/{{sample}}/repeatmasker/{{sample}}.fa.align",
        assembly_stats=f"{OUTDIR}/{{sample}}/genome/{{sample}}.assembly_stats.tsv",
    output:
        divsum=f"{OUTDIR}/{{arm}}/{{sample}}/divergence/{{sample}}.divsum",
        landscape=f"{OUTDIR}/{{arm}}/{{sample}}/divergence/{{sample}}.landscape.html",
    threads: config["resources"]["divergence"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["divergence"]["mem"] * attempt,
        hrs=config["resources"]["divergence"]["hrs"],
        shell_exec="bash",
    conda:
        "workflow/envs/repeatmasker.yaml"
    log:
        f"{OUTDIR}/logs/{{arm}}/{{sample}}/divergence.log",
    shell:
        """
        exec > {log} 2>&1
        {ENV_PATH_GUARD}
        mkdir -p $(dirname {output.divsum})

        # bioconda's RepeatMasker util scripts (calcDivergenceFromAlign.pl,
        # createRepeatLandscape.pl) ship with @INC missing the actual
        # RepeatMaskerConfig.pm directory (share/RepeatMasker itself, not
        # its util/ subdir) -- a known packaging gap, not a broken install.
        # Same defensive approach as setup_famdb's famdb.py search: find it
        # rather than guess a path that can vary by build/version.
        rm_config_path=$(find "$CONDA_PREFIX" -iname 'RepeatMaskerConfig.pm' 2>/dev/null | head -n1)
        if [ -z "$rm_config_path" ]; then
            echo "[ERROR] RepeatMaskerConfig.pm not found anywhere under \\$CONDA_PREFIX ($CONDA_PREFIX)." >&2
            echo "[ERROR] This conda RepeatMasker install may need (re)configuration -- see README's FamDB/RepeatMasker sections, or reinstall the env." >&2
            exit 1
        fi
        export PERL5LIB="$(dirname "$rm_config_path"):${{PERL5LIB:-}}"

        GENOME_BP=$(awk -F'\\t' 'NR==2{{print $3}}' {input.assembly_stats})
        calcDivergenceFromAlign.pl -s {output.divsum} -noCpGMod {input.align}
        createRepeatLandscape.pl -div {output.divsum} -g "$GENOME_BP" > {output.landscape}
        """


# -----------------------------------------------------------------------------
# 6.9 summarize (per arm, sample)
# -----------------------------------------------------------------------------
rule summarize:
    input:
        out_file=f"{OUTDIR}/{{arm}}/{{sample}}/repeatmasker/{{sample}}.fa.out",
        tbl_file=f"{OUTDIR}/{{arm}}/{{sample}}/repeatmasker/{{sample}}.fa.tbl",
        align_file=f"{OUTDIR}/{{arm}}/{{sample}}/repeatmasker/{{sample}}.fa.align",
        assembly_stats=f"{OUTDIR}/{{sample}}/genome/{{sample}}.assembly_stats.tsv",
        tandem_table=(
            f"{OUTDIR}/{{arm}}/{{sample}}/summary/family_tandem.tsv"
            if _as_bool(config["family_tandem"].get("carve_unknown", True))
            else []
        ),
        reclass_table=[f"{OUTDIR}/classify/unknown_reclassification.tsv"] if CLASSIFY_ON else [],
        curated_table=CURATED_TABLES,
    output:
        class_chunk=f"{OUTDIR}/{{arm}}/{{sample}}/summary/class_composition.tsv",
        family_chunk=f"{OUTDIR}/{{arm}}/{{sample}}/summary/family_composition.tsv",
        divergence_chunk=f"{OUTDIR}/{{arm}}/{{sample}}/summary/divergence_landscape.tsv",
    threads: config["resources"]["summarize"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["summarize"]["mem"] * attempt,
        hrs=config["resources"]["summarize"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/{{arm}}/{{sample}}/summarize.log",
    params:
        tissue=lambda wc: TISSUE_BY_SAMPLE[wc.sample],
        landscape_max_div=config["summary"]["landscape_max_div"],
        tandem_arg=lambda wc, input: f"--tandem-table {input.tandem_table}" if input.tandem_table else "",
        reclass_arg=lambda wc, input: f"--reclass-table {input.reclass_table}" if input.reclass_table else "",
        # user table first: the first table listing a family wins
        curated_arg=lambda wc, input: (
            f"--curated-table {input.curated_table}"
            + (f" --crosscheck-table {SATX_APPLIED}" if SATX_APPLY else "")
        ) if input.curated_table else "",
    shell:
        "python3 workflow/scripts/summarize_rm.py "
        "--out-file {input.out_file} --tbl-file {input.tbl_file} "
        "--align-file {input.align_file} --assembly-stats {input.assembly_stats} "
        "--arm {wildcards.arm} --sample {wildcards.sample} --tissue {params.tissue} "
        "--landscape-max-div {params.landscape_max_div} {params.tandem_arg} {params.reclass_arg} "
        "{params.curated_arg} --class-out {output.class_chunk} --family-out {output.family_chunk} "
        "--divergence-out {output.divergence_chunk} > {log} 2>&1"


# -----------------------------------------------------------------------------
# Unknown-family reclassification and class verification (classify.*): three
# cheap, independent screens of every consensus in the shared library -- TEsorter
# (protein domains), Rfam (structured RNA), DIAMOND against TE-free host
# proteins -- merged by reclassify_unknown.py into one table that summarize
# applies to .out/.align labels. Labels only: nothing is remasked.
# -----------------------------------------------------------------------------
rule extract_consensi:
    # Every family (the screens run once on the whole library: Unknown
    # families are reclassified, classified ones verified) + their classes.
    input:
        library=lambda wc: MASK_SHARED_LIBRARY or f"{OUTDIR}/library/shared_library.fa",
    output:
        fa=f"{OUTDIR}/classify/library_consensi.fa",
        classes=f"{OUTDIR}/classify/library_classes.tsv",
    threads: config["resources"]["classify_light"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["classify_light"]["mem"] * attempt,
        hrs=config["resources"]["classify_light"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/classify/extract_consensi.log",
    params:
        # De novo families only (prefixes added by prefix_library.py); curated
        # Dfam entries keep their labels unless classify.screen_dfam is set.
        keep=(
            ""
            if _as_bool(CLASSIFY.get("screen_dfam", False))
            else " ".join(f"--keep-prefix {sp}{config['library']['sample_prefix_sep']}" for sp in SAMPLE_IDS)
        ),
    shell:
        "python3 {SCRIPTS}/reclassify_unknown.py extract --library {input.library} --out {output.fa} "
        "--classes-out {output.classes} {params.keep} > {log} 2>&1"


rule tesorter_library:
    input:
        fa=f"{OUTDIR}/classify/library_consensi.fa",
    output:
        rexdb=f"{OUTDIR}/classify/tesorter/library.rexdb-metazoa.cls.tsv",
        gydb=f"{OUTDIR}/classify/tesorter/library.gydb.cls.tsv",
    threads: config["resources"]["tesorter"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["tesorter"]["mem"] * attempt,
        hrs=config["resources"]["tesorter"]["hrs"],
        shell_exec="bash",
    conda:
        "workflow/envs/tesorter.yaml"
    log:
        f"{OUTDIR}/logs/classify/tesorter_library.log",
    # -dp2: no TEsorter pass 2 (BLAST similarity to already-classified
    # sequences), so every call here is domain evidence. -nolib: we don't use
    # its RepeatMasker library.
    params:
        workdir=f"{OUTDIR}/classify/tesorter",
        fa_abs=lambda wc, input: os.path.abspath(input.fa),
        cov=CLASSIFY.get("tesorter_min_cov", 20),
        evalue=CLASSIFY.get("tesorter_max_evalue", 1e-3),
    shell:
        """
        exec > {log} 2>&1
        set -euo pipefail
        {ENV_PATH_GUARD}
        mkdir -p {params.workdir} && cd {params.workdir}
        for db in rexdb-metazoa gydb; do
            # TEsorter silently reuses a non-empty <pre>.domtbl from an earlier
            # run on a different FASTA (stale hits, or a KeyError on names no
            # longer present), so start each db clean.
            rm -rf library.$db.* tmp_$db
            TEsorter {params.fa_abs} -db $db -st nucl -p {threads} -cov {params.cov} -eval {params.evalue} \
                -dp2 -nolib -pre library.$db -tmp tmp_$db
            [ -f library.$db.cls.tsv ] || {{ echo "[ERROR] TEsorter wrote no library.$db.cls.tsv"; ls -l; exit 1; }}
            rm -rf tmp_$db
        done
        """


rule host_proteins_prep:
    input:
        annotation=ANNOT_PROTEINS,
        swissprot=[SWISSPROT] if SWISSPROT else [],
    output:
        annot=f"{OUTDIR}/classify/host/annotation_proteins.kwfiltered.faa",
        sprot=f"{OUTDIR}/classify/host/swissprot.kwfiltered.faa",
    threads: config["resources"]["classify_light"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["classify_light"]["mem"] * attempt,
        hrs=config["resources"]["classify_light"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/classify/host_proteins_prep.log",
    params:
        sprot_arg=f"--swissprot {SWISSPROT}" if SWISSPROT else "",
        keywords=" ".join(f"--keyword {shlex.quote(k)}" for k in CLASSIFY.get("te_protein_keywords", [])),
    shell:
        "python3 {SCRIPTS}/host_proteins.py prep --annotation {input.annotation} {params.sprot_arg} "
        "{params.keywords} --annot-out {output.annot} --sprot-out {output.sprot} > {log} 2>&1"


rule host_proteins_tesorter:
    # TE-domain screen of the annotation proteins (Swiss-Prot: keywords only).
    input:
        annot=f"{OUTDIR}/classify/host/annotation_proteins.kwfiltered.faa",
    output:
        cls=f"{OUTDIR}/classify/host/annotation.rexdb-metazoa.cls.tsv",
        dom=f"{OUTDIR}/classify/host/annotation.rexdb-metazoa.dom.tsv",
    threads: config["resources"]["tesorter"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["tesorter"]["mem"] * attempt,
        hrs=config["resources"]["tesorter"]["hrs"],
        shell_exec="bash",
    conda:
        "workflow/envs/tesorter.yaml"
    log:
        f"{OUTDIR}/logs/classify/host_proteins_tesorter.log",
    params:
        workdir=f"{OUTDIR}/classify/host",
        annot_abs=lambda wc, input: os.path.abspath(input.annot),
        cov=CLASSIFY.get("tesorter_min_cov", 20),
        evalue=CLASSIFY.get("tesorter_max_evalue", 1e-3),
    shell:
        """
        exec > {log} 2>&1
        set -euo pipefail
        {ENV_PATH_GUARD}
        cd {params.workdir}
        if [ ! -s {params.annot_abs} ]; then
            echo "no annotation proteins: nothing to screen"
            : > annotation.rexdb-metazoa.cls.tsv; : > annotation.rexdb-metazoa.dom.tsv
            exit 0
        fi
        rm -rf annotation.rexdb-metazoa.* tmp_prot  # no stale domtbl reuse (see tesorter_library)
        TEsorter {params.annot_abs} -db rexdb-metazoa -st prot -p {threads} -cov {params.cov} \
            -eval {params.evalue} -dp2 -nolib -pre annotation.rexdb-metazoa -tmp tmp_prot
        touch annotation.rexdb-metazoa.cls.tsv annotation.rexdb-metazoa.dom.tsv
        rm -rf tmp_prot
        """


rule host_protein_db:
    input:
        annot=f"{OUTDIR}/classify/host/annotation_proteins.kwfiltered.faa",
        sprot=f"{OUTDIR}/classify/host/swissprot.kwfiltered.faa",
        cls=f"{OUTDIR}/classify/host/annotation.rexdb-metazoa.cls.tsv",
        dom=f"{OUTDIR}/classify/host/annotation.rexdb-metazoa.dom.tsv",
    output:
        faa=f"{OUTDIR}/classify/host/host_proteins.faa",
        dmnd=f"{OUTDIR}/classify/host/host_proteins.dmnd",
    threads: config["resources"]["diamond"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["diamond"]["mem"] * attempt,
        hrs=config["resources"]["diamond"]["hrs"],
        shell_exec="bash",
    conda:
        "workflow/envs/diamond.yaml"
    log:
        f"{OUTDIR}/logs/classify/host_protein_db.log",
    params:
        db_prefix=lambda wc, output: output.dmnd[: -len(".dmnd")],
    shell:
        """
        exec > {log} 2>&1
        set -euo pipefail
        {ENV_PATH_GUARD}
        python3 {SCRIPTS}/host_proteins.py finalize --annot {input.annot} --sprot {input.sprot} \
            --tesorter {input.cls} {input.dom} --out {output.faa}
        diamond makedb --in {output.faa} -d {params.db_prefix} -p {threads}
        """


rule diamond_host:
    input:
        fa=f"{OUTDIR}/classify/library_consensi.fa",
        dmnd=f"{OUTDIR}/classify/host/host_proteins.dmnd",
    output:
        tsv=f"{OUTDIR}/classify/diamond_host.tsv",
    threads: config["resources"]["diamond"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["diamond"]["mem"] * attempt,
        hrs=config["resources"]["diamond"]["hrs"],
        shell_exec="bash",
    conda:
        "workflow/envs/diamond.yaml"
    log:
        f"{OUTDIR}/logs/classify/diamond_host.log",
    params:
        evalue=CLASSIFY.get("diamond_max_evalue", 1e-10),
    shell:
        """
        exec > {log} 2>&1
        set -euo pipefail
        {ENV_PATH_GUARD}
        diamond blastx -q {input.fa} -d {input.dmnd} -o {output.tsv} --more-sensitive \
            -e {params.evalue} -k 25 -p {threads} \
            --outfmt 6 qseqid sseqid pident length qstart qend qlen evalue bitscore stitle
        """


RFAM_PREP_CM = f"{OUTDIR}/classify/rfam/Rfam.cm"


rule rfam_prep:
    # Copy (or gunzip) the configured Rfam models into the run and cmpress
    # them in the pipeline's own Infernal env, so the user only downloads
    # Rfam.cm(.gz) -- nothing is written next to the original.
    input:
        cm=RFAM_CM if RFAM_CM else [],
    output:
        cm=RFAM_PREP_CM,
        idx=multiext(RFAM_PREP_CM, ".i1f", ".i1i", ".i1m", ".i1p"),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: config["resources"]["rfam_scan"]["mem"] * attempt,
        hrs=config["resources"]["rfam_scan"]["hrs"],
        shell_exec="bash",
    conda:
        "workflow/envs/infernal.yaml"
    log:
        f"{OUTDIR}/logs/classify/rfam_prep.log",
    shell:
        """
        exec > {log} 2>&1
        set -euo pipefail
        {ENV_PATH_GUARD}
        case {input.cm} in *.gz) gzip -dc {input.cm} > {output.cm} ;; *) cp {input.cm} {output.cm} ;; esac
        cmpress -F {output.cm}
        """


rule rfam_scan:
    input:
        fa=f"{OUTDIR}/classify/library_consensi.fa",
        cm=RFAM_PREP_CM,
        idx=multiext(RFAM_PREP_CM, ".i1f", ".i1i", ".i1m", ".i1p"),
    output:
        tblout=f"{OUTDIR}/classify/rfam.tblout",
    threads: config["resources"]["rfam_scan"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["rfam_scan"]["mem"] * attempt,
        hrs=config["resources"]["rfam_scan"]["hrs"],
        shell_exec="bash",
    conda:
        "workflow/envs/infernal.yaml"
    log:
        f"{OUTDIR}/logs/classify/rfam_scan.log",
    # No --clanin: it requires --fmt 2, which changes the tblout columns
    # reclassify_unknown.py reads, and the parser keeps only each consensus's
    # best-scoring hit anyway, so clan-overlap filtering adds nothing here.
    shell:
        """
        exec > {log} 2>&1
        set -euo pipefail
        {ENV_PATH_GUARD}
        cmscan --rfam --cut_ga --nohmmonly --cpu {threads} --tblout {output.tblout} \
            {input.cm} {input.fa} > /dev/null
        """


rule reclassify_unknown:
    input:
        fa=f"{OUTDIR}/classify/library_consensi.fa",
        classes=f"{OUTDIR}/classify/library_classes.tsv",
        rexdb=f"{OUTDIR}/classify/tesorter/library.rexdb-metazoa.cls.tsv",
        gydb=f"{OUTDIR}/classify/tesorter/library.gydb.cls.tsv",
        diamond=[f"{OUTDIR}/classify/diamond_host.tsv"] if HOST_SCREEN_ON else [],
        rfam=[f"{OUTDIR}/classify/rfam.tblout"] if RFAM_ON else [],
        # copy counts for classify.host_max_copies (TE ORFs annotated as genes)
        family_tandem=expand(f"{OUTDIR}/shared/{{sample}}/summary/family_tandem.tsv", sample=SAMPLE_IDS),
    output:
        tsv=f"{OUTDIR}/classify/unknown_reclassification.tsv",
    threads: config["resources"]["classify_light"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["classify_light"]["mem"] * attempt,
        hrs=config["resources"]["classify_light"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/classify/reclassify_unknown.log",
    params:
        diamond_arg=lambda wc, input: f"--diamond {input.diamond}" if input.diamond else "",
        rfam_arg=lambda wc, input: f"--rfam {input.rfam}" if input.rfam else "",
        min_domains=CLASSIFY.get("tesorter_min_domains", 1),
        min_host_cov=CLASSIFY.get("min_host_cov", 0.3),
        max_evalue=CLASSIFY.get("diamond_max_evalue", 1e-10),
        min_rfam_cov=CLASSIFY.get("min_rfam_cov", 0.5),
        host_max_copies=CLASSIFY.get("host_max_copies", 50),
        copy_scale=" ".join(f"{s}={c}" for s, c in COPIES_BY_SAMPLE.items() if c > 1),
    shell:
        "python3 {SCRIPTS}/reclassify_unknown.py merge --consensi {input.fa} --classes {input.classes} "
        "--tesorter-rexdb {input.rexdb} --tesorter-gydb {input.gydb} "
        "{params.diamond_arg} {params.rfam_arg} --min-domains {params.min_domains} "
        "--min-host-cov {params.min_host_cov} --max-evalue {params.max_evalue} "
        "--min-rfam-cov {params.min_rfam_cov} --family-tandem {input.family_tandem} "
        "--host-max-copies {params.host_max_copies} --copy-scale {params.copy_scale} --out {output.tsv} > {log} 2>&1"


DFAM_MATCHES = f"{OUTDIR}/library/dfam_overlap/dfam_matches.tsv"


def library_source_dfam_matches(wildcards):
    # summary/ (the full run) builds dfam_overlap first and tracks it.
    return [DFAM_MATCHES] if INCLUDE_DFAM and wildcards.sumdir == "summary" else []


def library_source_matches_arg(wildcards, input):
    # summary_shared_only/ reads the matches if they were already built but
    # never depends on them: an input edge (even ancient()) lets a stale
    # discovery chain be scheduled through dfam_overlap. Rebuilt matches
    # need --forcerun library_source.
    if input.dfam_matches:
        return f"--dfam-matches {input.dfam_matches}"
    if INCLUDE_DFAM and os.path.exists(DFAM_MATCHES):
        return f"--dfam-matches {DFAM_MATCHES}"
    return ""


rule library_source:
    # How much of each sample's masked bp comes from its own de novo
    # families, the other samples' de novo families, Dfam, and RepeatMasker's
    # built-in simple-repeat screen; de novo sources split into known-in-Dfam
    # vs novel when dfam_matches.tsv is available (README "Library sources").
    input:
        family_tandem=f"{OUTDIR}/{{sumdir}}/family_tandem.tsv",
        assembly_covariates=f"{OUTDIR}/{{sumdir}}/assembly_covariates.tsv",
        dfam_matches=library_source_dfam_matches,
    output:
        tsv=f"{OUTDIR}/{{sumdir}}/library_source.tsv",
    wildcard_constraints:
        sumdir="summary|summary_shared_only",
    threads: config["resources"]["classify_light"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["classify_light"]["mem"] * attempt,
        hrs=config["resources"]["classify_light"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/summary/library_source_{{sumdir}}.log",
    params:
        sample=" ".join(SAMPLE_IDS),
        sep=config["library"]["sample_prefix_sep"],
        matches_arg=library_source_matches_arg,
    shell:
        "python3 {SCRIPTS}/library_source.py --family-tandem {input.family_tandem} "
        "--assembly-covariates {input.assembly_covariates} --sample-ids {params.sample} "
        "--sep '{params.sep}' {params.matches_arg} --out {output.tsv} > {log} 2>&1"


rule element_groups:
    # Sum the bp of library families that are pieces of one element
    # (summary.element_groups, from family_groups.py members). README
    # "Element groups". With satellite_crosscheck.apply, the applied
    # proposals' groups (harmonized motif names) are reported too.
    input:
        groups=GROUP_TABLES + ([SATX_APPLIED] if SATX_APPLY else []),
        family_tandem=f"{OUTDIR}/{{sumdir}}/family_tandem.tsv",
        assembly_covariates=f"{OUTDIR}/{{sumdir}}/assembly_covariates.tsv",
    output:
        tsv=f"{OUTDIR}/{{sumdir}}/element_groups.tsv",
    wildcard_constraints:
        sumdir="summary|summary_shared_only",
    threads: 1
    resources:
        mem=lambda wildcards, attempt: config["resources"]["classify_light"]["mem"] * attempt,
        hrs=config["resources"]["classify_light"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/summary/element_groups_{{sumdir}}.log",
    params:
        sample=" ".join(SAMPLE_IDS),
    shell:
        "python3 {SCRIPTS}/family_groups.py report --groups {input.groups} "
        "--family-tandem {input.family_tandem} --assembly-covariates {input.assembly_covariates} "
        "--sample-ids {params.sample} --out {output.tsv} > {log} 2>&1"


rule curation_candidates:
    # Families worth a manual recheck: large and Unknown, consensus too long
    # for their label, flagged by verify_classes, or very young (README
    # "Curating a family").
    input:
        family_tandem=f"{OUTDIR}/{{sumdir}}/family_tandem.tsv",
        reclass=[f"{OUTDIR}/classify/unknown_reclassification.tsv"] if CLASSIFY_ON else [],
        disagreements=[f"{OUTDIR}/{{sumdir}}/class_disagreements.tsv"] if CLASSIFY_ON else [],
        groups=GROUP_TABLES + ([SATX_APPLIED] if SATX_APPLY else []),
    output:
        tsv=f"{OUTDIR}/{{sumdir}}/curation_candidates.tsv",
    wildcard_constraints:
        sumdir="summary|summary_shared_only",
    threads: 1
    resources:
        mem=lambda wildcards, attempt: config["resources"]["classify_light"]["mem"] * attempt,
        hrs=config["resources"]["classify_light"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/summary/curation_candidates_{{sumdir}}.log",
    params:
        sample=" ".join(SAMPLE_IDS),
        opt=lambda wc, input: " ".join(
            ([f"--reclass {input.reclass}"] if input.reclass else [])
            + ([f"--disagreements {input.disagreements}"] if input.disagreements else [])
            + ([f"--groups {input.groups}"] if input.groups else [])
        ),
        min_pct=config["summary"].get("curation_min_pct_masked", 0.1),
        young=config["summary"].get("curation_young_div", 3.0),
    shell:
        "python3 {SCRIPTS}/curation_candidates.py --family-tandem {input.family_tandem} {params.opt} "
        "--sample-ids {params.sample} --min-pct-masked {params.min_pct} --young-div {params.young} "
        "--out {output.tsv} > {log} 2>&1"


rule verify_classes:
    # Report-only check of RepeatClassifier's labels on the CLASSIFIED
    # families against the same screens, bp-weighted with the shared-arm
    # family_tandem tables. Nothing is relabelled (README).
    input:
        fa=f"{OUTDIR}/classify/library_consensi.fa",
        classes=f"{OUTDIR}/classify/library_classes.tsv",
        rexdb=f"{OUTDIR}/classify/tesorter/library.rexdb-metazoa.cls.tsv",
        gydb=f"{OUTDIR}/classify/tesorter/library.gydb.cls.tsv",
        diamond=[f"{OUTDIR}/classify/diamond_host.tsv"] if HOST_SCREEN_ON else [],
        rfam=[f"{OUTDIR}/classify/rfam.tblout"] if RFAM_ON else [],
        family_tandem=expand(f"{OUTDIR}/shared/{{sample}}/summary/family_tandem.tsv", sample=SAMPLE_IDS),
    output:
        verification=f"{OUTDIR}/{{sumdir}}/class_verification.tsv",
        disagreements=f"{OUTDIR}/{{sumdir}}/class_disagreements.tsv",
    wildcard_constraints:
        sumdir="summary|summary_shared_only",
    threads: config["resources"]["classify_light"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["classify_light"]["mem"] * attempt,
        hrs=config["resources"]["classify_light"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/classify/verify_classes_{{sumdir}}.log",
    params:
        diamond_arg=lambda wc, input: f"--diamond {input.diamond}" if input.diamond else "",
        rfam_arg=lambda wc, input: f"--rfam {input.rfam}" if input.rfam else "",
        min_domains=CLASSIFY.get("tesorter_min_domains", 1),
        min_host_cov=CLASSIFY.get("min_host_cov", 0.3),
        max_evalue=CLASSIFY.get("diamond_max_evalue", 1e-10),
        min_rfam_cov=CLASSIFY.get("min_rfam_cov", 0.5),
        host_max_copies=CLASSIFY.get("host_max_copies", 50),
        copy_scale=" ".join(f"{s}={c}" for s, c in COPIES_BY_SAMPLE.items() if c > 1),
    shell:
        "python3 {SCRIPTS}/reclassify_unknown.py verify --consensi {input.fa} --classes {input.classes} "
        "--tesorter-rexdb {input.rexdb} --tesorter-gydb {input.gydb} {params.diamond_arg} {params.rfam_arg} "
        "--family-tandem {input.family_tandem} --min-domains {params.min_domains} "
        "--min-host-cov {params.min_host_cov} --max-evalue {params.max_evalue} "
        "--min-rfam-cov {params.min_rfam_cov} --host-max-copies {params.host_max_copies} --copy-scale {params.copy_scale} "
        "--out {output.verification} "
        "--disagreements-out {output.disagreements} > {log} 2>&1"


# -----------------------------------------------------------------------------
# family_tandem (per arm, sample): is each library family arranged in
# tandem arrays (satellite-like) or dispersed? Same operational definition as
# the removed satellite arm's satellite_library_qc (README). summarize reads
# it to report tandem Unknown families as Unknown_tandem
# (family_tandem.carve_unknown). Combined tables are shared-arm only.
# -----------------------------------------------------------------------------
rule family_tandem:
    input:
        out_file=f"{OUTDIR}/{{arm}}/{{sample}}/repeatmasker/{{sample}}.fa.out",
        library=rm_library,
        assembly_stats=f"{OUTDIR}/{{sample}}/genome/{{sample}}.assembly_stats.tsv",
    output:
        family_chunk=f"{OUTDIR}/{{arm}}/{{sample}}/summary/family_tandem.tsv",
        class_chunk=f"{OUTDIR}/{{arm}}/{{sample}}/summary/class_tandem.tsv",
    threads: config["resources"]["family_tandem"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["family_tandem"]["mem"] * attempt,
        hrs=config["resources"]["family_tandem"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/{{arm}}/{{sample}}/family_tandem.log",
    params:
        t=config["family_tandem"],
        # dual_hap assemblies hold two copies of every locus: double the
        # absolute copy thresholds (percent-based thresholds are unaffected)
        min_copies=lambda wc: config["family_tandem"]["min_copies"] * COPIES_BY_SAMPLE[wc.sample],
        major_min_copies=lambda wc: config["family_tandem"]["major_min_copies"] * COPIES_BY_SAMPLE[wc.sample],
        major_min_bp=lambda wc: config["family_tandem"]["major_min_bp"] * COPIES_BY_SAMPLE[wc.sample],
    shell:
        "python3 workflow/scripts/family_tandem.py "
        "--out-file {input.out_file} --library {input.library} "
        "--assembly-stats {input.assembly_stats} "
        "--arm {wildcards.arm} --sample {wildcards.sample} "
        "--min-len {params.t[min_cons_len]} --max-len {params.t[max_cons_len]} "
        "--min-copies {params.min_copies} --min-array-copies {params.t[min_array_copies]} "
        "--min-tandem-frac {params.t[min_tandem_frac]} "
        "--max-short-period-frac {params.t[max_short_period_frac]} "
        "--major-min-copies {params.major_min_copies} --major-min-bp {params.major_min_bp} "
        "--out {output.family_chunk} --class-out {output.class_chunk} > {log} 2>&1"


rule combine_family_tandem:
    # Concatenate the per-sample tables (header once) into summary/ for a
    # full run or summary_shared_only/ for report_shared_only.
    input:
        family_chunks=expand(f"{OUTDIR}/shared/{{sample}}/summary/family_tandem.tsv", sample=SAMPLE_IDS),
        class_chunks=expand(f"{OUTDIR}/shared/{{sample}}/summary/class_tandem.tsv", sample=SAMPLE_IDS),
    output:
        family_tandem=f"{OUTDIR}/{{sumdir}}/family_tandem.tsv",
        class_tandem=f"{OUTDIR}/{{sumdir}}/class_tandem.tsv",
    wildcard_constraints:
        sumdir="summary|summary_shared_only",
    threads: 1
    resources:
        mem=lambda wildcards, attempt: config["resources"]["combine_summaries"]["mem"] * attempt,
        hrs=config["resources"]["combine_summaries"]["hrs"],
        shell_exec="bash",
    shell:
        "awk 'FNR == 1 && NR != 1 {{next}} {{print}}' {input.family_chunks} > {output.family_tandem} && "
        "awk 'FNR == 1 && NR != 1 {{next}} {{print}}' {input.class_chunks} > {output.class_tandem}"


rule round_saturation:
    # Own-arm masked bp by the RepeatModeler round that discovered each
    # family -- the data behind "should we sample more deeply?" (README).
    input:
        out_file=f"{OUTDIR}/own/{{sample}}/repeatmasker/{{sample}}.fa.out",
        assembly_stats=f"{OUTDIR}/{{sample}}/genome/{{sample}}.assembly_stats.tsv",
        tandem_table=f"{OUTDIR}/own/{{sample}}/summary/family_tandem.tsv",
    output:
        f"{OUTDIR}/own/{{sample}}/summary/round_saturation.tsv",
    threads: config["resources"]["summarize"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["summarize"]["mem"] * attempt,
        hrs=config["resources"]["summarize"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/own/{{sample}}/round_saturation.log",
    shell:
        "python3 {SCRIPTS}/round_saturation.py --out-file {input.out_file} "
        "--assembly-stats {input.assembly_stats} --sample {wildcards.sample} "
        "--tandem-table {input.tandem_table} --out {output} > {log} 2>&1"


rule round_novelty:
    # Per round: families that are new vs refinements of earlier families
    # (RepeatModeler's 80/80 cd-hit rule), with the own-arm bp each set owns.
    # The "another round?" number is new_dispersed_pct_non_n (README).
    input:
        rounds_fa=f"{OUTDIR}/{{sample}}/repeatmodeler/{{sample}}.rounds.consensi.fa",
        ltr_fa=f"{OUTDIR}/{{sample}}/ltr/{{sample}}.ltrs.fa",
        tandem_table=f"{OUTDIR}/own/{{sample}}/summary/family_tandem.tsv",
        assembly_stats=f"{OUTDIR}/{{sample}}/genome/{{sample}}.assembly_stats.tsv",
    output:
        f"{OUTDIR}/own/{{sample}}/summary/round_novelty.tsv",
    threads: config["resources"]["merge_families"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["merge_families"]["mem"] * attempt,
        hrs=config["resources"]["merge_families"]["hrs"],
        shell_exec="bash",
    singularity:
        TETOOLS
    log:
        f"{OUTDIR}/logs/own/{{sample}}/round_novelty.log",
    params:
        workdir=f"{OUTDIR}/own/{{sample}}/summary/round_novelty_work",
        rm_cfg=f"{SCRIPTS}/rm_config_path.sh",
    shell:
        """
        exec > {log} 2>&1
        set -euo pipefail
        CDHIT="$(bash {params.rm_cfg} CDHIT_DIR)/cd-hit-est"
        python3 {SCRIPTS}/round_novelty.py --rounds-fa {input.rounds_fa} --ltr-fa {input.ltr_fa} \
            --tandem-table {input.tandem_table} --assembly-stats {input.assembly_stats} \
            --sample {wildcards.sample} --cdhit "$CDHIT" --threads {threads} \
            --workdir {params.workdir} --out {output}
        rm -rf {params.workdir}
        """


rule combine_round_novelty:
    input:
        expand(f"{OUTDIR}/own/{{sample}}/summary/round_novelty.tsv", sample=SAMPLE_IDS),
    output:
        f"{OUTDIR}/summary/discovery_round_novelty.tsv",
    threads: 1
    resources:
        mem=lambda wildcards, attempt: config["resources"]["combine_summaries"]["mem"] * attempt,
        hrs=config["resources"]["combine_summaries"]["hrs"],
        shell_exec="bash",
    shell:
        "awk 'FNR == 1 && NR != 1 {{next}} {{print}}' {input} > {output}"


rule ltr_skipped_composition:
    # What the LTR tools' timed-out windows hold: shared-arm class and
    # tandem_family bp inside ltr/skipped_windows.tsv, vs genome-wide (README).
    input:
        skipped=f"{OUTDIR}/{{sample}}/ltr/skipped_windows.tsv",
        out_file=f"{OUTDIR}/shared/{{sample}}/repeatmasker/{{sample}}.fa.out",
        tandem_table=f"{OUTDIR}/shared/{{sample}}/summary/family_tandem.tsv",
        assembly_stats=f"{OUTDIR}/{{sample}}/genome/{{sample}}.assembly_stats.tsv",
    output:
        f"{OUTDIR}/shared/{{sample}}/summary/ltr_skipped_composition.tsv",
    threads: config["resources"]["summarize"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["summarize"]["mem"] * attempt,
        hrs=config["resources"]["summarize"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/shared/{{sample}}/ltr_skipped_composition.log",
    shell:
        "python3 {SCRIPTS}/ltr_skipped_composition.py --skipped {input.skipped} "
        "--out-file {input.out_file} --tandem-table {input.tandem_table} "
        "--assembly-stats {input.assembly_stats} --sample {wildcards.sample} "
        "--out {output} > {log} 2>&1"


# -----------------------------------------------------------------------------
# Satellite / rDNA / mito cross-check (satellite_crosscheck.*; README
# "Satellite, rDNA and mito cross-check"). Report-only: it reads existing
# shared-arm results and writes evidence + proposals; nothing is relabelled.
# Run on its own with the satellite_crosscheck_only target.
# -----------------------------------------------------------------------------
TC_CONSTRAINT = "|".join(re.escape(s) for s in TC_SAMPLES) or "__no_tidecluster_samples__"
SATDIR = f"{OUTDIR}/satellite_crosscheck"


rule trc_regions:
    # TideCluster arrays -> pipeline contig names, after checking TideCluster
    # ran on the same assembly (seqid_lengths.tsv vs fingerprint.tsv).
    input:
        name_map=f"{OUTDIR}/{{sample}}/genome/{{sample}}.name_map.tsv",
        fingerprint=f"{OUTDIR}/{{sample}}/genome/{{sample}}.fingerprint.tsv",
        clustering=lambda wc: _tc_file(wc.sample, "clustering.gff3"),
        seqid_lengths=lambda wc: _tc_file(wc.sample, "seqid_lengths.tsv"),
        consensus=lambda wc: _tc_file(wc.sample, "consensus_dimer_library.fasta"),
    output:
        regions=f"{OUTDIR}/shared/{{sample}}/satellite/trc_regions.tsv",
        trc_info=f"{OUTDIR}/shared/{{sample}}/satellite/trc_info.tsv",
        consensus=f"{OUTDIR}/shared/{{sample}}/satellite/trc_consensus.fa",
        params=f"{OUTDIR}/shared/{{sample}}/satellite/tidecluster_params.tsv",
    wildcard_constraints:
        sample=TC_CONSTRAINT,
    threads: config["resources"]["trc_regions"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["trc_regions"]["mem"] * attempt,
        hrs=config["resources"]["trc_regions"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/shared/{{sample}}/trc_regions.log",
    params:
        tc_dir=lambda wc: SATX_ANNOT[wc.sample]["tidecluster_dir"],
        tc_prefix=lambda wc: SATX_ANNOT[wc.sample]["tidecluster_prefix"],
        outdir=lambda wc, output: os.path.dirname(output.regions),
        allow=lambda wc: "--allow-seqid-mismatch" if _as_bool(SATX.get("allow_seqid_mismatch", False)) else "",
    shell:
        """
        exec > {log} 2>&1
        set -euo pipefail
        python3 {SCRIPTS}/tidecluster_regions.py --dir {params.tc_dir} --prefix {params.tc_prefix} \
            --sample {wildcards.sample} --name-map {input.name_map} --fingerprint {input.fingerprint} \
            {params.allow} --regions {output.regions} --trc-info {output.trc_info} \
            --consensus {output.consensus} --params {output.params}
        # TideCluster's own run record, for provenance
        for f in cmd_args.json pipeline_stats.json; do
            if [ -e {params.tc_dir}/{params.tc_prefix}_$f ]; then
                cp {params.tc_dir}/{params.tc_prefix}_$f {params.outdir}/tidecluster_$f
            fi
        done
        """


rule trc_crosscheck:
    # Which library families hold the TRC arrays, how much of each family
    # sits in them, and how much array sequence the library misses.
    input:
        out_file=f"{OUTDIR}/shared/{{sample}}/repeatmasker/{{sample}}.fa.out",
        tandem_table=f"{OUTDIR}/shared/{{sample}}/summary/family_tandem.tsv",
        regions=f"{OUTDIR}/shared/{{sample}}/satellite/trc_regions.tsv",
    output:
        trc=f"{OUTDIR}/shared/{{sample}}/satellite/trc_crosscheck.trc.tsv",
        family=f"{OUTDIR}/shared/{{sample}}/satellite/trc_crosscheck.family.tsv",
    wildcard_constraints:
        sample=TC_CONSTRAINT,
    threads: config["resources"]["trc_crosscheck"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["trc_crosscheck"]["mem"] * attempt,
        hrs=config["resources"]["trc_crosscheck"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/shared/{{sample}}/trc_crosscheck.log",
    shell:
        "python3 {SCRIPTS}/trc_crosscheck.py --out-file {input.out_file} "
        "--tandem-table {input.tandem_table} --regions {input.regions} --sample {wildcards.sample} "
        "--trc-out {output.trc} --family-out {output.family} > {log} 2>&1"


rule satellite_consensus_blast:
    # Candidate library consensi x TRC consensi x ribotin x mito, all-vs-all;
    # plus the whole library against ribotin + mito only.
    input:
        library=lambda wc: MASK_SHARED_LIBRARY or f"{OUTDIR}/library/shared_library.fa",
        family_tandem=f"{OUTDIR}/summary/family_tandem.tsv",
        crosscheck=expand(f"{OUTDIR}/shared/{{sample}}/satellite/trc_crosscheck.family.tsv", sample=TC_SAMPLES),
        trc_consensus=expand(f"{OUTDIR}/shared/{{sample}}/satellite/trc_consensus.fa", sample=TC_SAMPLES),
        refs=[r.split("=", 1)[1] for r in SATX_REFS],
        motifs=[SATX_MOTIFS] if SATX_MOTIFS else [],
    output:
        db_fa=f"{SATDIR}/blast/satellite_db.fa",
        refs_fa=f"{SATDIR}/blast/refs.fa",
        sat=f"{SATDIR}/blast/satellite_all_vs_all.tsv",
        ref=f"{SATDIR}/blast/library_vs_refs.tsv",
    threads: config["resources"]["satellite_blast"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["satellite_blast"]["mem"] * attempt,
        hrs=config["resources"]["satellite_blast"]["hrs"],
        shell_exec="bash",
    conda:
        "workflow/envs/blast.yaml"
    log:
        f"{OUTDIR}/logs/satellite_crosscheck/satellite_consensus_blast.log",
    params:
        refs=" ".join(SATX_REFS),
        focus=" ".join(SATX.get("focus_families") or []),
        min_trc_bp=SATX.get("candidate_min_trc_bp", 10000),
        evalue=SATX.get("blast_evalue", 1e-10),
        motifs=lambda wc, input: f"--motifs {input.motifs}" if input.motifs else "",
        max_targets=lambda wc, input: 2000 if input.motifs else 500,
        fmt="6 qseqid sseqid pident length qstart qend sstart send qlen slen evalue bitscore",
    shell:
        """
        exec > {log} 2>&1
        set -euo pipefail
        {ENV_PATH_GUARD}
        d=$(dirname {output.db_fa})
        python3 {SCRIPTS}/satellite_evidence.py select --library {input.library} \
            --family-tandem {input.family_tandem} --crosscheck {input.crosscheck} \
            --trc-consensus {input.trc_consensus} --ref {params.refs} --focus {params.focus} \
            --min-trc-bp {params.min_trc_bp} {params.motifs} \
            --db-fasta {output.db_fa} --refs-fasta {output.refs_fa}
        makeblastdb -in {output.db_fa} -dbtype nucl -out $d/satellite_db
        blastn -task blastn -query {output.db_fa} -db $d/satellite_db -evalue {params.evalue} \
            -dust no -num_threads {threads} -max_target_seqs {params.max_targets} -outfmt "{params.fmt}" -out {output.sat}
        if [ -s {output.refs_fa} ]; then
            makeblastdb -in {output.refs_fa} -dbtype nucl -out $d/refs
            blastn -task dc-megablast -query {input.library} -db $d/refs -evalue {params.evalue} \
                -num_threads {threads} -max_target_seqs 50 -outfmt "{params.fmt}" -out {output.ref}
        else
            : > {output.ref}
        fi
        """


rule satellite_evidence:
    input:
        family_tandem=f"{OUTDIR}/summary/family_tandem.tsv",
        crosscheck=expand(f"{OUTDIR}/shared/{{sample}}/satellite/trc_crosscheck.family.tsv", sample=TC_SAMPLES),
        trc_info=expand(f"{OUTDIR}/shared/{{sample}}/satellite/trc_info.tsv", sample=TC_SAMPLES),
        tc_params=expand(f"{OUTDIR}/shared/{{sample}}/satellite/tidecluster_params.tsv", sample=TC_SAMPLES),
        sat=f"{SATDIR}/blast/satellite_all_vs_all.tsv",
        ref=f"{SATDIR}/blast/library_vs_refs.tsv",
    output:
        evidence=f"{OUTDIR}/summary/satellite_family_evidence.tsv",
        calls=f"{OUTDIR}/summary/satellite_family_calls.tsv",
        pairs=f"{OUTDIR}/summary/satellite_family_pairs.tsv",
        proposals=f"{OUTDIR}/summary/satellite_proposals.tsv",
        trc_pairs=f"{OUTDIR}/summary/satellite_trc_pairs.tsv",
    threads: config["resources"]["satellite_evidence"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["satellite_evidence"]["mem"] * attempt,
        hrs=config["resources"]["satellite_evidence"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/satellite_crosscheck/satellite_evidence.log",
    params:
        infer=" ".join(f"{t}={d}" for t, d in SATX_INFER.items()),
        focus=" ".join(SATX.get("focus_families") or []),
        indep=" ".join(",".join(p) for p in (SATX.get("expected_independent") or [])),
        x=SATX,
        major_min_bp=config["family_tandem"]["major_min_bp"],
        copy_scale=" ".join(f"{s}={c}" for s, c in COPIES_BY_SAMPLE.items() if c > 1),
        eligible_min_share=SATX.get("eligible_min_share", 0.1),
        motif_args=" ".join(f"--{k.replace('_', '-')} {SATX[k]}" for k in
                            ("motif_min_id", "motif_high_cov", "motif_high_id", "motif_min_cov", "motif_partial_cov", "motif_array_tandem_frac")
                            if k in SATX),
    shell:
        "python3 {SCRIPTS}/satellite_evidence.py report --family-tandem {input.family_tandem} "
        "--crosscheck {input.crosscheck} --trc-info {input.trc_info} --tc-params {input.tc_params} "
        "--sat-blast {input.sat} --ref-blast {input.ref} --infer {params.infer} "
        "--focus {params.focus} --expected-independent {params.indep} "
        "--min-trc-cov {params.x[min_trc_cov]} --major-min-bp {params.major_min_bp} --copy-scale {params.copy_scale} "
        "--min-pair-id {params.x[min_pair_id]} --min-pair-cov {params.x[min_pair_cov]} "
        "--related-min-id {params.x[related_min_id]} --independent-max-cov {params.x[independent_max_cov]} "
        "--min-tandem-frac-long {params.x[min_tandem_frac_long]} --partial-trc-cov {params.x[partial_trc_cov]} "
        "--min-shared-trc-bp {params.x[min_shared_trc_bp]} "
        "--min-shared-trc-frac {params.x[min_shared_trc_frac]} "
        "--max-independent-shared-frac {params.x[max_independent_shared_frac]} "
        "--monomer-tol {params.x[monomer_tol]} "
        "--rdna-min-cov {params.x[rdna_min_cov]} --rdna-min-id {params.x[rdna_min_id]} "
        "--mito-min-cov {params.x[mito_min_cov]} --mito-min-id {params.x[mito_min_id]} "
        "--eligible-min-share {params.eligible_min_share} {params.motif_args} "
        "--evidence-out {output.evidence} --calls-out {output.calls} "
        "--pairs-out {output.pairs} --proposals-out {output.proposals} "
        "--trc-pairs-out {output.trc_pairs} > {log} 2>&1"


rule satellite_apply:
    # The cross-check proposals summarize applies (satellite_crosscheck.apply):
    # rows at apply_confidence, minus families the user's curated table lists.
    input:
        proposals=f"{OUTDIR}/summary/satellite_proposals.tsv",
        user=[CURATED_FAMILIES] if CURATED_FAMILIES else [],
    output:
        SATX_APPLIED,
    threads: config["resources"]["satellite_evidence"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["satellite_evidence"]["mem"] * attempt,
        hrs=config["resources"]["satellite_evidence"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/satellite_crosscheck/satellite_apply.log",
    params:
        user=lambda wc, input: f"--user-table {input.user}" if input.user else "",
        conf=" ".join(SATX.get("apply_confidence") or ["high"]),
    shell:
        "python3 {SCRIPTS}/satellite_evidence.py apply --proposals {input.proposals} {params.user} "
        "--confidence {params.conf} --out {output} > {log} 2>&1"


# Input-only target (no run/shell), like library_only: a target with a body
# becomes a cluster job and needs resources.
if "satellite_crosscheck_only" in sys.argv and not SATX_ON:
    print("[satellite_crosscheck] WARNING: satellite_crosscheck.external_annotations is not set; "
          "satellite_crosscheck_only has nothing to do.")


rule satellite_crosscheck_only:
    # Phase 1 entry point: the cross-check on existing results only. Run with
    # --rerun-triggers mtime so finished upstream jobs aren't redone.
    input:
        f"{OUTDIR}/summary/satellite_family_calls.tsv" if SATX_ON else [],
        expand(f"{OUTDIR}/shared/{{sample}}/satellite/trc_crosscheck.trc.tsv", sample=TC_SAMPLES),


# -----------------------------------------------------------------------------
# combine_summaries — final long-format comparison tables + arm_concordance
# (+ round saturation and LTR discovery summaries)
# -----------------------------------------------------------------------------
def _combine_inputs(wildcards):
    inputs = {
        "class_chunks": expand(
            f"{OUTDIR}/{{arm}}/{{sample}}/summary/class_composition.tsv", arm=ARMS, sample=SAMPLE_IDS
        ),
        "family_chunks": expand(
            f"{OUTDIR}/{{arm}}/{{sample}}/summary/family_composition.tsv", arm=ARMS, sample=SAMPLE_IDS
        ),
        "divergence_chunks": expand(
            f"{OUTDIR}/{{arm}}/{{sample}}/summary/divergence_landscape.tsv", arm=ARMS, sample=SAMPLE_IDS
        ),
        "assembly_stats_chunks": expand(
            f"{OUTDIR}/{{sample}}/genome/{{sample}}.assembly_stats.tsv", sample=SAMPLE_IDS
        ),
        "round_chunks": expand(f"{OUTDIR}/own/{{sample}}/summary/round_saturation.tsv", sample=SAMPLE_IDS),
        "ltr_summaries": expand(f"{OUTDIR}/{{sample}}/ltr/ltr_discovery_summary.tsv", sample=SAMPLE_IDS),
        "ltr_skipped": expand(f"{OUTDIR}/shared/{{sample}}/summary/ltr_skipped_composition.tsv", sample=SAMPLE_IDS),
    }
    return inputs


def _combine_outputs():
    outputs = {
        "class_composition": f"{OUTDIR}/summary/class_composition.tsv",
        "family_composition": f"{OUTDIR}/summary/family_composition.tsv",
        "divergence_landscape": f"{OUTDIR}/summary/divergence_landscape.tsv",
        "assembly_covariates": f"{OUTDIR}/summary/assembly_covariates.tsv",
        "arm_concordance": f"{OUTDIR}/summary/arm_concordance.tsv",
        "round_saturation": f"{OUTDIR}/summary/discovery_round_saturation.tsv",
        "ltr_discovery": f"{OUTDIR}/summary/ltr_discovery.tsv",
        "ltr_skipped": f"{OUTDIR}/summary/ltr_skipped_composition.tsv",
    }
    return outputs


rule combine_summaries:
    input:
        unpack(_combine_inputs),
    output:
        **_combine_outputs(),
    threads: config["resources"]["combine_summaries"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["combine_summaries"]["mem"] * attempt,
        hrs=config["resources"]["combine_summaries"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/summary/combine_summaries.log",
    shell:
        "python3 workflow/scripts/combine_summaries.py "
        "--manifest {config[manifest]} "
        "--class-chunks {input.class_chunks} "
        "--family-chunks {input.family_chunks} "
        "--divergence-chunks {input.divergence_chunks} "
        "--assembly-stats-chunks {input.assembly_stats_chunks} "
        "--class-composition-out {output.class_composition} "
        "--family-composition-out {output.family_composition} "
        "--divergence-landscape-out {output.divergence_landscape} "
        "--assembly-covariates-out {output.assembly_covariates} "
        "--arm-concordance-out {output.arm_concordance} "
        "--round-saturation-chunks {input.round_chunks} "
        "--round-saturation-out {output.round_saturation} "
        "--ltr-summary-chunks {input.ltr_summaries} "
        "--ltr-summary-out {output.ltr_discovery} "
        "--ltr-skipped-chunks {input.ltr_skipped} "
        "--ltr-skipped-out {output.ltr_skipped} "
        "> {log} 2>&1"


rule combine_summaries_shared_only:
    # combine_summaries for report_shared_only: shared-arm chunks only, and
    # none of the own-arm / discovery outputs.
    input:
        class_chunks=expand(f"{OUTDIR}/shared/{{sample}}/summary/class_composition.tsv", sample=SAMPLE_IDS),
        family_chunks=expand(f"{OUTDIR}/shared/{{sample}}/summary/family_composition.tsv", sample=SAMPLE_IDS),
        divergence_chunks=expand(f"{OUTDIR}/shared/{{sample}}/summary/divergence_landscape.tsv", sample=SAMPLE_IDS),
        assembly_stats_chunks=expand(f"{OUTDIR}/{{sample}}/genome/{{sample}}.assembly_stats.tsv", sample=SAMPLE_IDS),
    output:
        class_composition=f"{OUTDIR}/summary_shared_only/class_composition.tsv",
        family_composition=f"{OUTDIR}/summary_shared_only/family_composition.tsv",
        divergence_landscape=f"{OUTDIR}/summary_shared_only/divergence_landscape.tsv",
        assembly_covariates=f"{OUTDIR}/summary_shared_only/assembly_covariates.tsv",
    threads: config["resources"]["combine_summaries"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["combine_summaries"]["mem"] * attempt,
        hrs=config["resources"]["combine_summaries"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/summary_shared_only/combine_summaries.log",
    shell:
        "python3 workflow/scripts/combine_summaries.py "
        "--manifest {config[manifest]} "
        "--class-chunks {input.class_chunks} "
        "--family-chunks {input.family_chunks} "
        "--divergence-chunks {input.divergence_chunks} "
        "--assembly-stats-chunks {input.assembly_stats_chunks} "
        "--class-composition-out {output.class_composition} "
        "--family-composition-out {output.family_composition} "
        "--divergence-landscape-out {output.divergence_landscape} "
        "--assembly-covariates-out {output.assembly_covariates} "
        "> {log} 2>&1"


# -----------------------------------------------------------------------------
# provenance.txt
# -----------------------------------------------------------------------------
rule provenance:
    input:
        rm_version=f"{OUTDIR}/library/repeatmasker_version.txt",
        famdb_release=f"{OUTDIR}/library/famdb_release_info.txt",
        cdhit_version=f"{OUTDIR}/library/cdhit_version.txt",
        repeatmodeler_provenance=expand(
            f"{OUTDIR}/{{sample}}/repeatmodeler/{{sample}}.repeatmodeler_provenance.txt",
            sample=SAMPLE_IDS,
        ),
        fingerprints=expand(f"{OUTDIR}/{{sample}}/genome/{{sample}}.fingerprint.tsv", sample=SAMPLE_IDS),
        ltr_versions=expand(f"{OUTDIR}/{{sample}}/ltr/tool_versions.txt", sample=SAMPLE_IDS),
        finder_versions=expand(
            f"{OUTDIR}/{{sample}}/ltr/groups/{{group}}.finder.version.txt", sample=SAMPLE_IDS, group=LTR_GROUPS[:1]
        ) if USE_LTR_FINDER else [],
        ltr_discovery=f"{OUTDIR}/summary/ltr_discovery.tsv",
        config_snapshot="config.yaml",
    output:
        f"{OUTDIR}/summary/provenance.txt",
    threads: config["resources"]["provenance"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["provenance"]["mem"] * attempt,
        hrs=config["resources"]["provenance"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/summary/provenance.log",
    params:
        curated_override=config["library"]["curated_override"] or "(none — clustered library used)",
        container=TETOOLS,
    shell:
        """
        exec > {output} 2> {log}
        echo "=== Container (RepeatModeler, BuildDatabase, LTR pipeline, RepeatClassifier) ==="; echo "{params.container}"
        echo; echo "=== RepeatMasker ==="; cat {input.rm_version}
        echo; echo "=== FamDB / Dfam release ==="; cat {input.famdb_release}
        echo; echo "=== cd-hit ==="; cat {input.cdhit_version}
        echo; echo "=== RepeatModeler version + in-container Dfam partitions, per sample ==="
        for f in {input.repeatmodeler_provenance}; do echo "--- $f ---"; cat "$f"; done
        echo; echo "=== Genome fingerprints (md5 / n_seqs / total_bp) ==="
        for f in {input.fingerprints}; do echo "--- $f ---"; grep '^#' "$f"; done
        echo; echo "=== LTR discovery tool versions, per sample ==="
        for f in {input.ltr_versions} {input.finder_versions}; do echo "--- $f ---"; cat "$f"; done
        echo "vendored: see workflow/vendor/README.md (LTR_HARVEST_parallel c3c9b3c, LTR_FINDER_parallel f1036ca, RepeatModeler 2.0.9 LTRPipeline, all patched)"
        echo; echo "=== LTR candidates and window-timeout skipped bp ==="; cat {input.ltr_discovery}
        echo; echo "=== curated_override ==="; echo "{params.curated_override}"
        echo; echo "=== config.yaml snapshot ==="; cat {input.config_snapshot}
        """


# -----------------------------------------------------------------------------
# 6.10 plot
# -----------------------------------------------------------------------------
rule plot:
    input:
        class_composition=f"{OUTDIR}/summary/class_composition.tsv",
        divergence_landscape=f"{OUTDIR}/summary/divergence_landscape.tsv",
        arm_concordance=f"{OUTDIR}/summary/arm_concordance.tsv",
        assembly_covariates=f"{OUTDIR}/summary/assembly_covariates.tsv",
    output:
        class_plot=f"{OUTDIR}/plots/class_composition_shared.png",
        divergence_plot=f"{OUTDIR}/plots/divergence_landscape.png",
        concordance_plot=f"{OUTDIR}/plots/arm_concordance.png",
    threads: config["resources"]["plot"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["plot"]["mem"] * attempt,
        hrs=config["resources"]["plot"]["hrs"],
        shell_exec="bash",
    conda:
        "workflow/envs/r_plot.yaml"
    params:
        out_dir=f"{OUTDIR}/plots",
        sample_order=PLOT_SAMPLE_ORDER,
    log:
        f"{OUTDIR}/logs/summary/plot.log",
    shell:
        "Rscript workflow/scripts/plot_repeat_compare.R "
        "{input.class_composition} {input.divergence_landscape} "
        "{input.arm_concordance} {input.assembly_covariates} {params.out_dir} "
        "{params.sample_order} > {log} 2>&1"


rule plot_shared_only:
    # plot for report_shared_only: no arm_concordance (no own arm).
    input:
        class_composition=f"{OUTDIR}/summary_shared_only/class_composition.tsv",
        divergence_landscape=f"{OUTDIR}/summary_shared_only/divergence_landscape.tsv",
        assembly_covariates=f"{OUTDIR}/summary_shared_only/assembly_covariates.tsv",
    output:
        class_plot=f"{OUTDIR}/plots_shared_only/class_composition_shared.png",
        divergence_plot=f"{OUTDIR}/plots_shared_only/divergence_landscape.png",
    threads: config["resources"]["plot"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["plot"]["mem"] * attempt,
        hrs=config["resources"]["plot"]["hrs"],
        shell_exec="bash",
    conda:
        "workflow/envs/r_plot.yaml"
    params:
        out_dir=f"{OUTDIR}/plots_shared_only",
        sample_order=PLOT_SAMPLE_ORDER,
    log:
        f"{OUTDIR}/logs/summary_shared_only/plot.log",
    shell:
        "Rscript workflow/scripts/plot_repeat_compare.R "
        "{input.class_composition} {input.divergence_landscape} "
        "NONE {input.assembly_covariates} {params.out_dir} {params.sample_order} "
        "> {log} 2>&1"
