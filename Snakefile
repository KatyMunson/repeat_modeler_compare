# =============================================================================
# Snakefile — repeat_compare: RepeatModeler2 + RepeatMasker cross-species
# repeat comparison
#
# 1. Sanitizes each manifest species' genome FASTA and computes assembly QC
#    covariates (contig N50, N content, etc — the denominators used
#    everywhere downstream).
# 2. Runs RepeatModeler2 (via the dfam/tetools Singularity image) de novo on
#    each species independently.
# 3. Prefixes each species' family names with its species_id (RepeatModeler
#    names collide across species), optionally exports a Dfam supplement,
#    and clusters all species' families together with cd-hit-est into one
#    non-redundant "shared" library — the primary, apples-to-apples masking
#    library. A per-species "own" library (no cross-species clustering) is
#    also assembled as a secondary sanity-check arm.
# 4. Masks every species' genome with BOTH the shared library and its own
#    library (two "arms"), computes divergence landscapes, and summarizes
#    non-overlapping repeat bp per class and per Class/Family, reported
#    against both the total and non-N assembly length.
# 5. Combines everything into long-format comparison tables and plots.
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
                raise ValueError(
                    f"{path}:{lineno}: expected 5 tab-separated fields "
                    f"(species_id, species_name, fasta, tissue, accession), "
                    f"got {len(fields)}: {line!r}"
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

OUTDIR = config["outdir"]
ARMS = ["shared", "own"]
INCLUDE_DFAM = bool(config["library"]["include_dfam"])
DFAM_TAXON = config["library"]["dfam_taxon"]
DFAM_EXPORT_FASTA = f"{OUTDIR}/library/dfam_{DFAM_TAXON}.fa"
LTRSTRUCT_FLAG = "-LTRStruct" if config["repeatmodeler"]["ltrstruct"] else ""

wildcard_constraints:
    species="|".join(re.escape(s) for s in SPECIES_IDS),
    arm="shared|own",


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
        config["repeatmodeler"]["container"]
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
# 6.4 repeatmodeler
# -----------------------------------------------------------------------------
rule repeatmodeler:
    input:
        db_done=f"{OUTDIR}/{{species}}/repeatmodeler/{{species}}.build_db.done",
    output:
        families=f"{OUTDIR}/{{species}}/repeatmodeler/{{species}}-families.fa",
        stk=f"{OUTDIR}/{{species}}/repeatmodeler/{{species}}-families.stk",
        provenance=f"{OUTDIR}/{{species}}/repeatmodeler/{{species}}.repeatmodeler_provenance.txt",
    threads: config["resources"]["repeatmodeler"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["repeatmodeler"]["mem"] * attempt,
        hrs=config["resources"]["repeatmodeler"]["hrs"],
        shell_exec="bash",
    singularity:
        config["repeatmodeler"]["container"]
    params:
        workdir=f"{OUTDIR}/{{species}}/repeatmodeler",
        ltrstruct=LTRSTRUCT_FLAG,
        extra_args=config["repeatmodeler"]["extra_args"],
    log:
        f"{OUTDIR}/logs/{{species}}/repeatmodeler.log",
    shell:
        """
        exec > {log} 2>&1
        cd {params.workdir}

        # Restart safety (§6.4): reuse a previous RM_* working dir on retry
        # rather than starting over, instead of temp()-wiping it.
        RM_DIR=$(ls -d RM_* 2>/dev/null | head -n1 || true)
        if [ -n "$RM_DIR" ]; then
            RECOVER="-recoverDir $RM_DIR"
            echo "[INFO] Recovering previous RepeatModeler run from $RM_DIR"
        else
            RECOVER=""
        fi

        RepeatModeler -database {wildcards.species} -threads {threads} \
            {params.ltrstruct} $RECOVER {params.extra_args}

        {{
            echo "=== RepeatModeler version ==="
            RepeatModeler -version 2>&1 || echo "RepeatModeler -version failed"
            echo
            echo "=== famdb.py info (as seen by RepeatClassifier inside this container) ==="
            famdb.py info 2>&1 || echo "famdb.py info failed or not found in container"
        }} > $(basename {output.provenance})
        """


# -----------------------------------------------------------------------------
# 6.5 prefix_library
# -----------------------------------------------------------------------------
rule prefix_library:
    input:
        fa=f"{OUTDIR}/{{species}}/repeatmodeler/{{species}}-families.fa",
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
            "famdb.py families -f fasta_name --include-class-in-name -a -d "
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
    shell:
        "python3 workflow/scripts/library_membership.py "
        "--clstr {input.clstr} --sep {params.sep} --out {output} > {log} 2>&1"


def _shared_library_inputs(wildcards):
    inputs = {"nr": f"{OUTDIR}/library/shared_denovo.nr.fa"}
    if INCLUDE_DFAM:
        inputs["dfam"] = DFAM_EXPORT_FASTA
    return inputs


rule assemble_shared_library:
    # Primary comparison-arm library. If library.curated_override is set,
    # it replaces clustering entirely for THIS file (§8) — cluster_library
    # and library_membership.tsv above still run unconditionally, since the
    # de novo shared-vocabulary comparison is a result in its own right
    # independent of whether a curated library is used for masking.
    input:
        unpack(_shared_library_inputs),
    output:
        f"{OUTDIR}/library/shared_library.fa",
    threads: config["resources"]["assemble_shared_library"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["assemble_shared_library"]["mem"] * attempt,
        hrs=config["resources"]["assemble_shared_library"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/library/assemble_shared_library.log",
    params:
        curated_override=config["library"]["curated_override"],
    run:
        import shutil

        with open(log[0], "w") as logf:
            if params.curated_override:
                logf.write(f"Using curated_override: {params.curated_override}\n")
                shutil.copy(params.curated_override, output[0])
            else:
                logf.write("Building shared_library.fa from cluster_library + dfam export\n")
                with open(output[0], "w") as out_f:
                    with open(input.nr) as in_f:
                        shutil.copyfileobj(in_f, out_f)
                    if INCLUDE_DFAM:
                        with open(input.dfam) as in_f:
                            shutil.copyfileobj(in_f, out_f)


def _own_library_inputs(wildcards):
    inputs = {"prefixed": f"{OUTDIR}/{wildcards.species}/library/{wildcards.species}.prefixed.fa"}
    if INCLUDE_DFAM:
        inputs["dfam"] = DFAM_EXPORT_FASTA
    return inputs


rule own_library:
    # Sanity-check arm library: that species' de novo families only (+ same
    # optional Dfam export), assembled in its own small rule so both arms
    # call RepeatMasker identically (§6.7).
    input:
        unpack(_own_library_inputs),
    output:
        f"{OUTDIR}/own/{{species}}/library/{{species}}.own_library.fa",
    threads: config["resources"]["own_library"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["own_library"]["mem"] * attempt,
        hrs=config["resources"]["own_library"]["hrs"],
        shell_exec="bash",
    log:
        f"{OUTDIR}/logs/{{species}}/own_library.log",
    run:
        import shutil

        with open(output[0], "w") as out_f:
            with open(input.prefixed) as in_f:
                shutil.copyfileobj(in_f, out_f)
            if INCLUDE_DFAM:
                with open(input.dfam) as in_f:
                    shutil.copyfileobj(in_f, out_f)


# -----------------------------------------------------------------------------
# 6.8 repeatmasker (both arms)
# -----------------------------------------------------------------------------
def rm_library(wildcards):
    if wildcards.arm == "shared":
        return f"{OUTDIR}/library/shared_library.fa"
    return f"{OUTDIR}/own/{wildcards.species}/library/{wildcards.species}.own_library.fa"


rule repeatmasker:
    input:
        fa=f"{OUTDIR}/{{species}}/genome/{{species}}.fa",
        lib=rm_library,
        famdb_verified=f"{OUTDIR}/library/famdb_verified.txt",
    output:
        out_file=f"{OUTDIR}/{{arm}}/{{species}}/repeatmasker/{{species}}.fa.out",
        tbl_file=f"{OUTDIR}/{{arm}}/{{species}}/repeatmasker/{{species}}.fa.tbl",
        align_file=f"{OUTDIR}/{{arm}}/{{species}}/repeatmasker/{{species}}.fa.align",
    threads: config["resources"]["repeatmasker"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["repeatmasker"]["mem"] * attempt,
        hrs=config["resources"]["repeatmasker"]["hrs"],
        shell_exec="bash",
    conda:
        "workflow/envs/repeatmasker.yaml"
    params:
        outdir=f"{OUTDIR}/{{arm}}/{{species}}/repeatmasker",
        pa=config["resources"]["repeatmasker"]["threads"] // config["repeatmasker"]["cores_per_pa"],
        sensitive_flag="-s" if config["repeatmasker"]["sensitive"] else "",
        extra_args=config["repeatmasker"]["extra_args"],
        # Resolved here (DAG-build time, original invocation directory) --
        # NOT via a shell-level `readlink -f` inside the rule body, which
        # would run after the `cd {params.outdir}` below and resolve these
        # relative paths against the wrong directory, silently expanding to
        # nothing (readlink -f fails when the leading path components don't
        # exist, and a failed command substitution doesn't trip `set -e`).
        fa_abs=lambda wc, input: os.path.abspath(input.fa),
        lib_abs=lambda wc, input: os.path.abspath(input.lib),
    log:
        f"{OUTDIR}/logs/{{arm}}/{{species}}/repeatmasker.log",
    shell:
        "mkdir -p {params.outdir} && "
        "(cd {params.outdir} && trap 'rm -rf RM_*' EXIT && "
        "RepeatMasker -pa {params.pa} -lib {params.lib_abs} -xsmall -gff -a "
        "{params.sensitive_flag} {params.extra_args} -dir . {params.fa_abs}) "
        "> {log} 2>&1 && "
        "touch {output.out_file} {output.tbl_file} {output.align_file}"


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


# -----------------------------------------------------------------------------
# combine_summaries — final long-format comparison tables + arm_concordance
# -----------------------------------------------------------------------------
rule combine_summaries:
    input:
        class_chunks=expand(
            f"{OUTDIR}/{{arm}}/{{species}}/summary/class_composition.tsv", arm=ARMS, species=SPECIES_IDS
        ),
        family_chunks=expand(
            f"{OUTDIR}/{{arm}}/{{species}}/summary/family_composition.tsv", arm=ARMS, species=SPECIES_IDS
        ),
        divergence_chunks=expand(
            f"{OUTDIR}/{{arm}}/{{species}}/summary/divergence_landscape.tsv", arm=ARMS, species=SPECIES_IDS
        ),
        assembly_stats_chunks=expand(
            f"{OUTDIR}/{{species}}/genome/{{species}}.assembly_stats.tsv", species=SPECIES_IDS
        ),
    output:
        class_composition=f"{OUTDIR}/summary/class_composition.tsv",
        family_composition=f"{OUTDIR}/summary/family_composition.tsv",
        divergence_landscape=f"{OUTDIR}/summary/divergence_landscape.tsv",
        assembly_covariates=f"{OUTDIR}/summary/assembly_covariates.tsv",
        arm_concordance=f"{OUTDIR}/summary/arm_concordance.tsv",
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
    shell:
        """
        exec > {output} 2> {log}
        echo "=== RepeatMasker ==="; cat {input.rm_version}
        echo; echo "=== FamDB / Dfam release ==="; cat {input.famdb_release}
        echo; echo "=== cd-hit ==="; cat {input.cdhit_version}
        echo; echo "=== RepeatModeler version + in-container Dfam partitions, per species ==="
        for f in {input.repeatmodeler_provenance}; do echo "--- $f ---"; cat "$f"; done
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
