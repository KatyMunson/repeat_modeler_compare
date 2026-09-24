# =============================================================================
# Snakefile — repeat_compare: RepeatModeler2 + RepeatMasker cross-species
# repeat comparison
#
# 1. Sanitizes each manifest species' genome FASTA, computes assembly QC
#    covariates (contig N50, N content, etc — the denominators used
#    everywhere downstream) and a genome fingerprint (identity guard).
# 2. Runs RepeatModeler2's RECON/RepeatScout rounds (dfam/tetools
#    Singularity image) de novo on each species independently — WITHOUT
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
# 4. Prefixes each species' family names with its species_id, optionally
#    exports a Dfam supplement, clusters all species' families together
#    with cd-hit-est into one non-redundant "shared" library (+ Dfam), and
#    assembles a per-species "own" library.
# 5. Masks every species' genome with BOTH libraries (two "arms"), computes
#    divergence landscapes, and summarizes non-overlapping repeat bp per
#    class and per Class/Family against both total and non-N length.
# 6. Combines everything into long-format comparison tables and plots.
#
# A satellite screen / satellite-library arm existed briefly and was removed
# pending fixes to the satellite caller (compare_assemblies_satellites); it
# is preserved at commit 502accf (tag satellite-arm-v1) -- see README.
#
# Species undergoing programmed germline-to-soma genome rearrangement (e.g.
# hagfish) can have very different repeat content between tissues — this
# pipeline does not act on the manifest's `tissue` column, it only carries
# it into every summary table and warns if species disagree or are
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

configfile: "config.yaml"


# -----------------------------------------------------------------------------
# Manifest parsing (stdlib only, no pandas) — §5.1
# -----------------------------------------------------------------------------
def parse_manifest(path):
    manifest = []
    seen_ids = set()
    with open(path) as fh:
        for lineno, raw_line in enumerate(fh, 1):
            line = raw_line.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 5:
                extra = (" (a 6th satellite_lib column is no longer supported -- the satellite "
                         "arm was removed; see README)") if len(fields) == 6 else ""
                raise ValueError(
                    f"{path}:{lineno}: expected 5 tab-separated fields "
                    f"(species_id, species_name, fasta, tissue, accession), "
                    f"got {len(fields)}{extra}: {line!r}"
                )
            species_id, species_name, fasta, tissue, accession = fields
            if not re.fullmatch(r"[A-Za-z0-9]+", species_id):
                raise ValueError(
                    f"{path}:{lineno}: species_id '{species_id}' must match [A-Za-z0-9]+"
                )
            if species_id in seen_ids:
                raise ValueError(f"{path}:{lineno}: duplicate species_id '{species_id}'")
            seen_ids.add(species_id)
            if tissue not in ("germline", "soma", "unknown"):
                raise ValueError(
                    f"{path}:{lineno}: tissue must be one of germline|soma|unknown, got '{tissue}'"
                )
            if not os.path.exists(fasta):
                raise ValueError(f"{path}:{lineno}: fasta path does not exist: {fasta}")
            manifest.append(
                {
                    "species_id": species_id,
                    "species_name": species_name,
                    "fasta": fasta,
                    "tissue": tissue,
                    "accession": accession,
                }
            )
    if len(manifest) < 2:
        raise ValueError(f"{path}: at least 2 species are required, found {len(manifest)}")
    return manifest


MANIFEST = parse_manifest(config["manifest"])
SPECIES_IDS = [row["species_id"] for row in MANIFEST]
FASTA_BY_SPECIES = {row["species_id"]: row["fasta"] for row in MANIFEST}
TISSUE_BY_SPECIES = {row["species_id"]: row["tissue"] for row in MANIFEST}

tissue_values = set(TISSUE_BY_SPECIES.values())
if "unknown" in tissue_values or len(tissue_values) > 1:
    print(
        f"[repeat_compare] WARNING: manifest tissue values are {sorted(tissue_values)} — "
        f"a germline assembly and a somatic assembly are different genomes; this is not "
        f"treated as an error, but check before comparing across species.",
    )

# repeatmasker.sensitive is a single config value used by every arm/species by
# construction; validate its type here so a future refactor can't silently
# make it per-arm/per-species without this assertion catching it (§6.8).
if not isinstance(config["repeatmasker"]["sensitive"], bool):
    raise ValueError("config['repeatmasker']['sensitive'] must be true or false")

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
SCRIPTS = "workflow/scripts"
VENDOR = "workflow/vendor"


LTR_CFG = config["ltr_discovery"]
LTR_GROUPS = [f"g{i}" for i in range(int(LTR_CFG["n_groups"]))]
if config["resources"]["ltr_harvest_group"]["threads"] < 2 or config["resources"]["ltr_finder_group"]["threads"] < 2:
    # The vendored *_parallel scripts' 1-thread branch ignores -size/-time
    # (whole group, no timeout) -- exactly the stall this side pipeline avoids.
    raise ValueError("resources.ltr_harvest_group/ltr_finder_group threads must be >= 2")
USE_LTR_FINDER = bool(LTR_CFG["use_ltr_finder"])
LTR_TOOLS = ["harvest", "finder"] if USE_LTR_FINDER else ["harvest"]

wildcard_constraints:
    species="|".join(re.escape(s) for s in SPECIES_IDS),
    arm="shared|own",
    group=r"g\d+",
    tool="harvest|finder",

# Chunk count for the repeatmasker scatter/gather split (split_genome /
# repeatmasker_chunk / gather_repeatmasker below) -- reused pattern from
# compare_assemblies_satellites' stage 03 (its -lib-mode RepeatMasker
# scatter/gather), adapted here since our own repeatmasker rule originally
# ran each species' whole genome as one unchunked job.
scattergather:
    genome_chunks=config["repeatmasker"]["scatter_count"],


# -----------------------------------------------------------------------------
# rule all / library_only — §6.11
# -----------------------------------------------------------------------------
rule all:
    input:
        expand(
            f"{OUTDIR}/{{arm}}/{{species}}/repeatmasker/{{species}}.fa.out",
            arm=ARMS,
            species=SPECIES_IDS,
        ),
        f"{OUTDIR}/summary/class_composition.tsv",
        f"{OUTDIR}/summary/family_composition.tsv",
        f"{OUTDIR}/summary/divergence_landscape.tsv",
        f"{OUTDIR}/summary/arm_concordance.tsv",
        f"{OUTDIR}/summary/assembly_covariates.tsv",
        f"{OUTDIR}/summary/discovery_round_saturation.tsv",
        f"{OUTDIR}/summary/ltr_discovery.tsv",
        f"{OUTDIR}/summary/provenance.txt",
        f"{OUTDIR}/library/library_membership.tsv",
        f"{OUTDIR}/plots/class_composition_shared.png",
        f"{OUTDIR}/plots/divergence_landscape.png",
        f"{OUTDIR}/plots/arm_concordance.png",


rule library_only:
    # Stops after shared_library.fa so it can be inspected / manually
    # curated (§8) before the expensive masking step.
    input:
        f"{OUTDIR}/library/shared_library.fa",
        f"{OUTDIR}/library/library_membership.tsv",


# -----------------------------------------------------------------------------
# 6.1 prep_genome
# -----------------------------------------------------------------------------
rule prep_genome:
    input:
        fasta=lambda wc: FASTA_BY_SPECIES[wc.species],
    output:
        fa=f"{OUTDIR}/{{species}}/genome/{{species}}.fa",
        name_map=f"{OUTDIR}/{{species}}/genome/{{species}}.name_map.tsv",
    threads: config["resources"]["prep_genome"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["prep_genome"]["mem"] * attempt,
        hrs=config["resources"]["prep_genome"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/{{species}}/prep_genome.log",
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
        fa=f"{OUTDIR}/{{species}}/genome/{{species}}.fa",
    output:
        f"{OUTDIR}/{{species}}/genome/{{species}}.fingerprint.tsv",
    threads: config["resources"]["genome_fingerprint"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["genome_fingerprint"]["mem"] * attempt,
        hrs=config["resources"]["genome_fingerprint"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/{{species}}/genome_fingerprint.log",
    shell:
        "python3 {SCRIPTS}/fingerprint.py write --fasta {input.fa} --out {output} > {log} 2>&1"


# -----------------------------------------------------------------------------
# 6.2 assembly_stats
# -----------------------------------------------------------------------------
rule assembly_stats:
    input:
        fa=f"{OUTDIR}/{{species}}/genome/{{species}}.fa",
    output:
        f"{OUTDIR}/{{species}}/genome/{{species}}.assembly_stats.tsv",
    threads: config["resources"]["assembly_stats"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["assembly_stats"]["mem"] * attempt,
        hrs=config["resources"]["assembly_stats"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/{{species}}/assembly_stats.log",
    shell:
        "python3 workflow/scripts/assembly_stats.py "
        "--fasta {input.fa} --species-id {wildcards.species} --out {output} "
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
        fa=f"{OUTDIR}/{{species}}/genome/{{species}}.fa",
    output:
        touch(f"{OUTDIR}/{{species}}/repeatmodeler/{{species}}.build_db.done"),
    threads: config["resources"]["build_db"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["build_db"]["mem"] * attempt,
        hrs=config["resources"]["build_db"]["hrs"],
        shell_exec="bash",
    singularity:
        TETOOLS
    params:
        workdir=f"{OUTDIR}/{{species}}/repeatmodeler",
    log:
        f"{OUTDIR}/logs/{{species}}/build_db.log",
    shell:
        "mkdir -p {params.workdir} && "
        # This container's RepeatModeler 2.0.9 BuildDatabase rejects -engine
        # outright ("Unknown option: engine") -- WU-BLAST support was
        # dropped upstream and NCBI/RMBlast is the only engine now, so the
        # flag is no longer accepted at all, not just unnecessary.
        "BuildDatabase -name {params.workdir}/{wildcards.species} {input.fa} "
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
#  - A finished directory (consensi.fa.classified present) is used as-is.
#  - Otherwise -recoverDir. RepeatModeler's recovery exits 0 WITHOUT doing
#    anything ("appears to contain a successful run") when every round
#    already completed; since only the rounds are needed here, that case
#    is accepted too.
# -----------------------------------------------------------------------------
rule repeatmodeler:
    input:
        db_done=f"{OUTDIR}/{{species}}/repeatmodeler/{{species}}.build_db.done",
        fingerprint=f"{OUTDIR}/{{species}}/genome/{{species}}.fingerprint.tsv",
    output:
        consensi=f"{OUTDIR}/{{species}}/repeatmodeler/{{species}}.rounds.consensi.fa",
        stk=f"{OUTDIR}/{{species}}/repeatmodeler/{{species}}.rounds.families.stk",
        provenance=f"{OUTDIR}/{{species}}/repeatmodeler/{{species}}.repeatmodeler_provenance.txt",
    threads: config["resources"]["repeatmodeler"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["repeatmodeler"]["mem"] * attempt,
        hrs=config["resources"]["repeatmodeler"]["hrs"],
        shell_exec="bash",
    singularity:
        TETOOLS
    params:
        workdir=f"{OUTDIR}/{{species}}/repeatmodeler",
        extra_args=config["repeatmodeler"]["extra_args"],
        fp_abs=lambda wc, input: os.path.abspath(input.fingerprint),
        script_abs=os.path.abspath(f"{SCRIPTS}/fingerprint.py"),
        out_consensi=lambda wc, output: os.path.abspath(output.consensi),
        out_stk=lambda wc, output: os.path.abspath(output.stk),
        out_prov=lambda wc, output: os.path.abspath(output.provenance),
    log:
        f"{OUTDIR}/logs/{{species}}/repeatmodeler.log",
    shell:
        """
        exec > {log} 2>&1
        set -euo pipefail
        cd {params.workdir}

        RM_DIR=$(ls -d RM_* 2>/dev/null | head -n1 || true)
        if [ -z "$RM_DIR" ]; then
            cp {params.fp_abs} rm_run.fingerprint.tsv
            RepeatModeler -database {wildcards.species} -threads {threads} {params.extra_args}
            RM_DIR=$(ls -d RM_* | head -n1)
        else
            if [ ! -s rm_run.fingerprint.tsv ]; then
                echo "[ERROR] $RM_DIR exists but rm_run.fingerprint.tsv doesn't -- it wasn't started by this rule."
                echo "[ERROR] Move $RM_DIR out of $(pwd) (or delete it) and rerun."
                exit 1
            fi
            python3 {params.script_abs} check --expected rm_run.fingerprint.tsv --observed {params.fp_abs} || {{
                echo "[ERROR] $RM_DIR was started on a different genome than the current prepped FASTA."
                echo "[ERROR] Refusing to -recoverDir it. Move $RM_DIR and rm_run.fingerprint.tsv away and rerun."
                exit 1
            }}
            if [ -s "$RM_DIR/consensi.fa.classified" ]; then
                echo "[INFO] $RM_DIR already finished; reusing its rounds output"
            else
                echo "[INFO] Recovering previous RepeatModeler run from $RM_DIR"
                RepeatModeler -database {wildcards.species} -threads {threads} \
                    -recoverDir "$RM_DIR" {params.extra_args} | tee recover.stdout
                if [ ! -s "$RM_DIR/consensi.fa.classified" ] && \
                   ! grep -q "appears to contain a successful run" recover.stdout; then
                    echo "[ERROR] RepeatModeler recovery did not complete; see above."
                    exit 1
                fi
            fi
        fi

        test -s "$RM_DIR/consensi.fa" && test -s "$RM_DIR/families.stk"
        cp "$RM_DIR/consensi.fa" {params.out_consensi}
        cp "$RM_DIR/families.stk" {params.out_stk}
        {{
            echo "=== RepeatModeler version ==="
            RepeatModeler -version 2>&1 || echo "RepeatModeler -version failed"
            echo "=== run: rounds only (no -LTRStruct), extra_args='{params.extra_args}', dir $RM_DIR ==="
            echo
            echo "=== famdb.py info (as seen by RepeatClassifier inside this container) ==="
            famdb.py info 2>&1 || echo "famdb.py info failed or not found in container"
        }} > {params.out_prov}
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
        f"{OUTDIR}/{{species}}/genome/{{species}}.fa",
    output:
        groups=temp(expand(f"{OUTDIR}/{{{{species}}}}/ltr/groups/{{group}}.fa", group=LTR_GROUPS)),
        manifest=f"{OUTDIR}/{{species}}/ltr/groups/manifest.tsv",
    threads: config["resources"]["ltr_group_genome"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["ltr_group_genome"]["mem"] * attempt,
        hrs=config["resources"]["ltr_group_genome"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/{{species}}/ltr_group_genome.log",
    shell:
        "python3 {SCRIPTS}/group_genome.py --fasta {input} --outputs {output.groups} "
        "--manifest {output.manifest} > {log} 2>&1"


_EMPTY_SCN_HEADER = "# no sequences in this group"


rule ltr_harvest_group:
    input:
        f"{OUTDIR}/{{species}}/ltr/groups/{{group}}.fa",
    output:
        scn=f"{OUTDIR}/{{species}}/ltr/groups/{{group}}.harvest.scn",
        timeouts=f"{OUTDIR}/{{species}}/ltr/groups/{{group}}.harvest.timeouts.tsv",
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
        f"{OUTDIR}/logs/{{species}}/ltr_harvest_group/{{group}}.log",
    params:
        workdir=f"{OUTDIR}/{{species}}/ltr/work/harvest_{{group}}",
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
        f"{OUTDIR}/{{species}}/ltr/groups/{{group}}.fa",
    output:
        scn=f"{OUTDIR}/{{species}}/ltr/groups/{{group}}.finder.scn",
        timeouts=f"{OUTDIR}/{{species}}/ltr/groups/{{group}}.finder.timeouts.tsv",
        version=f"{OUTDIR}/{{species}}/ltr/groups/{{group}}.finder.version.txt",
    threads: config["resources"]["ltr_finder_group"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["ltr_finder_group"]["mem"] * attempt,
        hrs=config["resources"]["ltr_finder_group"]["hrs"],
        shell_exec="bash",
    retries: 1
    conda:
        "workflow/envs/ltr_finder.yaml"
    log:
        f"{OUTDIR}/logs/{{species}}/ltr_finder_group/{{group}}.log",
    params:
        workdir=f"{OUTDIR}/{{species}}/ltr/work/finder_{{group}}",
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
        genome=f"{OUTDIR}/{{species}}/genome/{{species}}.fa",
        harvest=expand(f"{OUTDIR}/{{{{species}}}}/ltr/groups/{{group}}.harvest.scn", group=LTR_GROUPS),
        harvest_to=expand(f"{OUTDIR}/{{{{species}}}}/ltr/groups/{{group}}.harvest.timeouts.tsv", group=LTR_GROUPS),
        finder=expand(f"{OUTDIR}/{{{{species}}}}/ltr/groups/{{group}}.finder.scn", group=LTR_GROUPS) if USE_LTR_FINDER else [],
        finder_to=expand(f"{OUTDIR}/{{{{species}}}}/ltr/groups/{{group}}.finder.timeouts.tsv", group=LTR_GROUPS) if USE_LTR_FINDER else [],
    output:
        scn=f"{OUTDIR}/{{species}}/ltr/rawLTR.scn",
        skipped=f"{OUTDIR}/{{species}}/ltr/skipped_windows.tsv",
        summary=f"{OUTDIR}/{{species}}/ltr/ltr_discovery_summary.tsv",
    threads: config["resources"]["ltr_gather"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["ltr_gather"]["mem"] * attempt,
        hrs=config["resources"]["ltr_gather"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/{{species}}/ltr_gather.log",
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
        "--window-size {params.size} --overlap {params.overlap} --species {wildcards.species} "
        "--out-scn {output.scn} --out-skipped {output.skipped} --out-summary {output.summary} "
        "> {log} 2>&1"


rule ltr_pipeline:
    input:
        genome=f"{OUTDIR}/{{species}}/genome/{{species}}.fa",
        scn=f"{OUTDIR}/{{species}}/ltr/rawLTR.scn",
    output:
        fa=f"{OUTDIR}/{{species}}/ltr/{{species}}.ltrs.fa",
        stk=f"{OUTDIR}/{{species}}/ltr/{{species}}.ltrs.stk",
        versions=f"{OUTDIR}/{{species}}/ltr/tool_versions.txt",
    threads: config["resources"]["ltr_pipeline"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["ltr_pipeline"]["mem"] * attempt,
        hrs=config["resources"]["ltr_pipeline"]["hrs"],
        shell_exec="bash",
    singularity:
        TETOOLS
    log:
        f"{OUTDIR}/logs/{{species}}/ltr_pipeline.log",
    params:
        workdir=f"{OUTDIR}/{{species}}/ltr/work/pipeline",
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
        ln -s {params.genome_abs} {wildcards.species}.fa
        export RM_DIR=$(dirname "$(readlink -f "$(command -v RepeatModeler)")")
        perl {params.tool} -inscn {params.scn_abs} -threads {threads} -tmpdir . {wildcards.species}.fa
        if [ -s {wildcards.species}.fa-ltrs.fa ]; then
            cp {wildcards.species}.fa-ltrs.fa {params.fa_abs}
            cp {wildcards.species}.fa-ltrs.stk {params.stk_abs}
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
        rounds_fa=f"{OUTDIR}/{{species}}/repeatmodeler/{{species}}.rounds.consensi.fa",
        rounds_stk=f"{OUTDIR}/{{species}}/repeatmodeler/{{species}}.rounds.families.stk",
        ltr_fa=f"{OUTDIR}/{{species}}/ltr/{{species}}.ltrs.fa",
        ltr_stk=f"{OUTDIR}/{{species}}/ltr/{{species}}.ltrs.stk",
    output:
        fa=f"{OUTDIR}/{{species}}/families/{{species}}.merged.consensi.fa",
        stk=f"{OUTDIR}/{{species}}/families/{{species}}.merged.families.stk",
    threads: config["resources"]["merge_families"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["merge_families"]["mem"] * attempt,
        hrs=config["resources"]["merge_families"]["hrs"],
        shell_exec="bash",
    singularity:
        TETOOLS
    log:
        f"{OUTDIR}/logs/{{species}}/merge_families.log",
    params:
        workdir=f"{OUTDIR}/{{species}}/families/merge_work",
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
        fa=f"{OUTDIR}/{{species}}/families/{{species}}.merged.consensi.fa",
        stk=f"{OUTDIR}/{{species}}/families/{{species}}.merged.families.stk",
    output:
        fa=f"{OUTDIR}/{{species}}/families/{{species}}-families.fa",
        stk=f"{OUTDIR}/{{species}}/families/{{species}}-families.stk",
    threads: config["resources"]["classify_families"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["classify_families"]["mem"] * attempt,
        hrs=config["resources"]["classify_families"]["hrs"],
        shell_exec="bash",
    singularity:
        TETOOLS
    log:
        f"{OUTDIR}/logs/{{species}}/classify_families.log",
    params:
        workdir=f"{OUTDIR}/{{species}}/families/classify_work",
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
        fa=f"{OUTDIR}/{{species}}/families/{{species}}-families.fa",
    output:
        f"{OUTDIR}/{{species}}/library/{{species}}.prefixed.fa",
    threads: config["resources"]["prefix_library"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["prefix_library"]["mem"] * attempt,
        hrs=config["resources"]["prefix_library"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/{{species}}/prefix_library.log",
    params:
        sep=config["library"]["species_prefix_sep"],
    shell:
        "python3 workflow/scripts/prefix_library.py "
        "--fasta {input.fa} --species-id {wildcards.species} --sep {params.sep} "
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
            # the curated cross-species reference set that section actually
            # describes.
            "famdb.py families -f fasta_name --include-class-in-name -c -a -d "
            "--add-reverse-complement '{params.taxon}' > {output} 2> {log}"


# -----------------------------------------------------------------------------
# 6.7 cluster_library + library_membership + shared/own library assembly
# -----------------------------------------------------------------------------
rule cluster_library:
    input:
        prefixed=expand(f"{OUTDIR}/{{species}}/library/{{species}}.prefixed.fa", species=SPECIES_IDS),
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
        sep=config["library"]["species_prefix_sep"],
        codes=" ".join(SPECIES_IDS),
        dfam_arg=lambda wc, input: f"--dfam {input.dfam}" if input.dfam else "",
    shell:
        "python3 workflow/scripts/library_membership.py "
        "--clstr {input.clstr} --sep {params.sep} --species-codes {params.codes} "
        "{params.dfam_arg} --out {output} > {log} 2>&1"


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
    inputs = {"prefixed": f"{OUTDIR}/{wildcards.species}/library/{wildcards.species}.prefixed.fa"}
    if INCLUDE_DFAM:
        inputs["dfam"] = DFAM_EXPORT_FASTA
    return inputs


rule own_library:
    # Sanity-check arm library: that species' de novo families only (+ the
    # same optional Dfam export), assembled in its own small rule so both
    # arms call RepeatMasker identically (§6.7).
    input:
        unpack(_own_library_inputs),
    output:
        fa=f"{OUTDIR}/own/{{species}}/library/{{species}}.own_library.fa",
        report=f"{OUTDIR}/own/{{species}}/library/{{species}}.own_library.sources.tsv",
    threads: config["resources"]["own_library"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["own_library"]["mem"] * attempt,
        hrs=config["resources"]["own_library"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/{{species}}/own_library.log",
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
# gather_repeatmasker_genome_lib), rather than running each species' whole
# genome as one unchunked job: real multi-Gb assemblies schedule much better
# on SGE as N independently-restartable chunk jobs. Deliberately NOT copying
# that stage's -nolow flag -- it suppresses RepeatMasker's built-in
# low-complexity/simple-repeat screen, which we need (Simple_repeat and
# Low_complexity are both required canonical output classes here).
# -----------------------------------------------------------------------------
def rm_library(wildcards):
    if wildcards.arm == "shared":
        return f"{OUTDIR}/library/shared_library.fa"
    return f"{OUTDIR}/own/{wildcards.species}/library/{wildcards.species}.own_library.fa"


rule split_genome:
    # Scatters one species' genome into config["repeatmasker"]["scatter_count"]
    # chunks (workflow/scripts/split_fasta.py, copied verbatim from the
    # sibling repo's common/scripts/split_fasta.py -- snake/boustrophedon by
    # contig count). Per-species, not per-arm: the genome being split is
    # identical regardless of which library later masks it.
    input:
        fasta=f"{OUTDIR}/{{species}}/genome/{{species}}.fa",
    output:
        fasta=temp(scatter.genome_chunks(
            f"{OUTDIR}/{{{{species}}}}/genome/chunks/{{scatteritem}}/{{scatteritem}}.fa"
        )),
    threads: config["resources"]["split_genome"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["split_genome"]["mem"] * attempt,
        hrs=config["resources"]["split_genome"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/{{species}}/split_genome.log",
    shell:
        "python3 workflow/scripts/split_fasta.py --infile {input.fasta} "
        "--outputs {output.fasta} > {log} 2>&1"


rule repeatmasker_chunk:
    input:
        fasta=f"{OUTDIR}/{{species}}/genome/chunks/{{scatteritem}}/{{scatteritem}}.fa",
        lib=rm_library,
        famdb_verified=f"{OUTDIR}/library/famdb_verified.txt",
    output:
        out_file=temp(f"{OUTDIR}/{{arm}}/{{species}}/repeatmasker/chunks/{{scatteritem}}/{{scatteritem}}.fa.out"),
        tbl_file=temp(f"{OUTDIR}/{{arm}}/{{species}}/repeatmasker/chunks/{{scatteritem}}/{{scatteritem}}.fa.tbl"),
        align_file=temp(f"{OUTDIR}/{{arm}}/{{species}}/repeatmasker/chunks/{{scatteritem}}/{{scatteritem}}.fa.align"),
    threads: config["resources"]["repeatmasker"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["repeatmasker"]["mem"] * attempt,
        hrs=config["resources"]["repeatmasker"]["hrs"],
        shell_exec="bash",
    conda:
        "workflow/envs/repeatmasker.yaml"
    params:
        outdir=f"{OUTDIR}/{{arm}}/{{species}}/repeatmasker/chunks/{{scatteritem}}",
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
        f"{OUTDIR}/logs/{{arm}}/{{species}}/repeatmasker_chunk/{{scatteritem}}.log",
    shell:
        "mkdir -p {params.outdir} && "
        "(cd {params.outdir} && trap 'rm -rf RM_*' EXIT && "
        "RepeatMasker -pa {params.pa} -lib {params.lib_abs} -xsmall -gff -a "
        "{params.sensitive_flag} {params.extra_args} -dir . {params.fa_abs}) "
        "> {log} 2>&1 && "
        "touch {output.out_file} {output.tbl_file} {output.align_file}"


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
            f"{OUTDIR}/{{{{arm}}}}/{{{{species}}}}/repeatmasker/chunks/{{scatteritem}}/{{scatteritem}}.fa.out"
        ),
        tbl_chunks=gather.genome_chunks(
            f"{OUTDIR}/{{{{arm}}}}/{{{{species}}}}/repeatmasker/chunks/{{scatteritem}}/{{scatteritem}}.fa.tbl"
        ),
        align_chunks=gather.genome_chunks(
            f"{OUTDIR}/{{{{arm}}}}/{{{{species}}}}/repeatmasker/chunks/{{scatteritem}}/{{scatteritem}}.fa.align"
        ),
    output:
        out_file=f"{OUTDIR}/{{arm}}/{{species}}/repeatmasker/{{species}}.fa.out",
        tbl_file=f"{OUTDIR}/{{arm}}/{{species}}/repeatmasker/{{species}}.fa.tbl",
        align_file=f"{OUTDIR}/{{arm}}/{{species}}/repeatmasker/{{species}}.fa.align",
    threads: config["resources"]["gather_repeatmasker"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["gather_repeatmasker"]["mem"] * attempt,
        hrs=config["resources"]["gather_repeatmasker"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/{{arm}}/{{species}}/gather_repeatmasker.log",
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
# -----------------------------------------------------------------------------
rule divergence:
    input:
        align=f"{OUTDIR}/{{arm}}/{{species}}/repeatmasker/{{species}}.fa.align",
        assembly_stats=f"{OUTDIR}/{{species}}/genome/{{species}}.assembly_stats.tsv",
    output:
        divsum=f"{OUTDIR}/{{arm}}/{{species}}/divergence/{{species}}.divsum",
        landscape=f"{OUTDIR}/{{arm}}/{{species}}/divergence/{{species}}.landscape.html",
    threads: config["resources"]["divergence"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["divergence"]["mem"] * attempt,
        hrs=config["resources"]["divergence"]["hrs"],
        shell_exec="bash",
    conda:
        "workflow/envs/repeatmasker.yaml"
    log:
        f"{OUTDIR}/logs/{{arm}}/{{species}}/divergence.log",
    shell:
        """
        exec > {log} 2>&1
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
# 6.9 summarize (per arm, species)
# -----------------------------------------------------------------------------
rule summarize:
    input:
        out_file=f"{OUTDIR}/{{arm}}/{{species}}/repeatmasker/{{species}}.fa.out",
        tbl_file=f"{OUTDIR}/{{arm}}/{{species}}/repeatmasker/{{species}}.fa.tbl",
        divsum=f"{OUTDIR}/{{arm}}/{{species}}/divergence/{{species}}.divsum",
        assembly_stats=f"{OUTDIR}/{{species}}/genome/{{species}}.assembly_stats.tsv",
    output:
        class_chunk=f"{OUTDIR}/{{arm}}/{{species}}/summary/class_composition.tsv",
        family_chunk=f"{OUTDIR}/{{arm}}/{{species}}/summary/family_composition.tsv",
        divergence_chunk=f"{OUTDIR}/{{arm}}/{{species}}/summary/divergence_landscape.tsv",
    threads: config["resources"]["summarize"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["summarize"]["mem"] * attempt,
        hrs=config["resources"]["summarize"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/{{arm}}/{{species}}/summarize.log",
    params:
        tissue=lambda wc: TISSUE_BY_SPECIES[wc.species],
        landscape_max_div=config["summary"]["landscape_max_div"],
    shell:
        "python3 workflow/scripts/summarize_rm.py "
        "--out-file {input.out_file} --tbl-file {input.tbl_file} "
        "--divsum-file {input.divsum} --assembly-stats {input.assembly_stats} "
        "--arm {wildcards.arm} --species {wildcards.species} --tissue {params.tissue} "
        "--landscape-max-div {params.landscape_max_div} "
        "--class-out {output.class_chunk} --family-out {output.family_chunk} "
        "--divergence-out {output.divergence_chunk} > {log} 2>&1"


rule round_saturation:
    # Own-arm masked bp by the RepeatModeler round that discovered each
    # family -- the data behind "should we sample more deeply?" (README).
    input:
        out_file=f"{OUTDIR}/own/{{species}}/repeatmasker/{{species}}.fa.out",
        assembly_stats=f"{OUTDIR}/{{species}}/genome/{{species}}.assembly_stats.tsv",
    output:
        f"{OUTDIR}/own/{{species}}/summary/round_saturation.tsv",
    threads: config["resources"]["summarize"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["summarize"]["mem"] * attempt,
        hrs=config["resources"]["summarize"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/own/{{species}}/round_saturation.log",
    shell:
        "python3 {SCRIPTS}/round_saturation.py --out-file {input.out_file} "
        "--assembly-stats {input.assembly_stats} --species {wildcards.species} "
        "--out {output} > {log} 2>&1"


# -----------------------------------------------------------------------------
# combine_summaries — final long-format comparison tables + arm_concordance
# (+ round saturation and LTR discovery summaries)
# -----------------------------------------------------------------------------
def _combine_inputs(wildcards):
    inputs = {
        "class_chunks": expand(
            f"{OUTDIR}/{{arm}}/{{species}}/summary/class_composition.tsv", arm=ARMS, species=SPECIES_IDS
        ),
        "family_chunks": expand(
            f"{OUTDIR}/{{arm}}/{{species}}/summary/family_composition.tsv", arm=ARMS, species=SPECIES_IDS
        ),
        "divergence_chunks": expand(
            f"{OUTDIR}/{{arm}}/{{species}}/summary/divergence_landscape.tsv", arm=ARMS, species=SPECIES_IDS
        ),
        "assembly_stats_chunks": expand(
            f"{OUTDIR}/{{species}}/genome/{{species}}.assembly_stats.tsv", species=SPECIES_IDS
        ),
        "round_chunks": expand(f"{OUTDIR}/own/{{species}}/summary/round_saturation.tsv", species=SPECIES_IDS),
        "ltr_summaries": expand(f"{OUTDIR}/{{species}}/ltr/ltr_discovery_summary.tsv", species=SPECIES_IDS),
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
            f"{OUTDIR}/{{species}}/repeatmodeler/{{species}}.repeatmodeler_provenance.txt",
            species=SPECIES_IDS,
        ),
        fingerprints=expand(f"{OUTDIR}/{{species}}/genome/{{species}}.fingerprint.tsv", species=SPECIES_IDS),
        ltr_versions=expand(f"{OUTDIR}/{{species}}/ltr/tool_versions.txt", species=SPECIES_IDS),
        finder_versions=expand(
            f"{OUTDIR}/{{species}}/ltr/groups/{{group}}.finder.version.txt", species=SPECIES_IDS, group=LTR_GROUPS[:1]
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
        echo; echo "=== RepeatModeler version + in-container Dfam partitions, per species ==="
        for f in {input.repeatmodeler_provenance}; do echo "--- $f ---"; cat "$f"; done
        echo; echo "=== Genome fingerprints (md5 / n_seqs / total_bp) ==="
        for f in {input.fingerprints}; do echo "--- $f ---"; grep '^#' "$f"; done
        echo; echo "=== LTR discovery tool versions, per species ==="
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
    log:
        f"{OUTDIR}/logs/summary/plot.log",
    shell:
        "Rscript workflow/scripts/plot_repeat_compare.R "
        "{input.class_composition} {input.divergence_landscape} "
        "{input.arm_concordance} {input.assembly_covariates} {params.out_dir} "
        "> {log} 2>&1"
