# repeat_compare

RepeatModeler2 + RepeatMasker cross-species repeat comparison, built for
*Eptatretus stoutii* (de novo) vs *Myxine limosa* (published), but
species-agnostic and manifest-driven.

## Design

Both genomes are masked with **two arms**:

| arm | library | role |
|---|---|---|
| `shared` | non-redundant union of all species' de novo RepeatModeler2 families (+ optional Dfam export) | **primary comparison** |
| `own` | that species' de novo families only (+ same optional Dfam export) | sanity check / concordance |

Masking every species with only its own library biases the comparison
(each genome is best-annotated for its own families), so `shared` is the
one to trust for cross-species numbers; `own` exists to sanity-check it via
`results/summary/arm_concordance.tsv`.

Hagfish undergo programmed germline-to-soma genome rearrangement — a
germline assembly and a somatic assembly are different genomes. This
pipeline does **not** act on the manifest's `tissue` column; it only
carries it into every summary table and warns (never fails) if species
disagree or are `unknown`. Assembly-quality covariates (contig count,
total/N/non-N length, N50) are reported alongside every repeat result for
the same reason — a fragmented or collapsed assembly undercounts repeats.

## Quickstart

```bash
# 1. Edit manifest.tsv with real absolute FASTA paths and confirmed tissue values.
# 2. Dry run (catches DAG/wildcard/syntax errors, not real module/tool availability):
snakemake -s Snakefile --configfile config.yaml -n -p --restart-times 3

# 3. Wiring test on a small subsample (set genome_prep.test_subsample_bp in config.yaml,
#    e.g. to 5000000, first) via a real cluster submission:
./runsnake 8 --configfile config.yaml all

# 4. Library-only target, to inspect/curate the shared library before the
#    expensive masking step:
./runsnake 8 --configfile config.yaml library_only
```

**`--configfile`-before-targets caveat:** passing a target name immediately
after `--configfile` makes Snakemake treat it as a second configfile.
Always put targets *before* `--configfile`, or use a `--` separator —
`./runsnake 8 --configfile config.yaml all` (configfile before the `all`
target) is correct; `./runsnake 8 all --configfile config.yaml` is not.

## Cluster gotcha: `snakemake/9.3.0` module vs. `--use-conda`

On liger, the `snakemake/9.3.0` modulefile **hard-conflicts with any
`miniconda` module** (`module show snakemake/9.3.0` lists `conflict
miniconda`) and, separately, unconditionally prepends its own bundled
conda (4.12.0, i.e. `conda 22.9.0`) onto `$PATH` itself — confirmed via
`which conda` after `module load snakemake/9.3.0` with no miniconda module
loaded at all. That conda is too old for Snakemake 9.x's `--use-conda`,
which requires **conda ≥24.7.1**, and since the module's own PATH prepend
happens at load time, it wins over any newer conda you put on `$PATH`
*before* `./runsnake` runs `module load snakemake/9.3.0` internally.
`module swap` doesn't help either, since there's no miniconda module
actually loaded to swap out.

The fix: install a personal, newer conda (or mamba/miniforge) anywhere in
your home directory, then point `runsnake` at its `bin/` directory via the
`CONDA_OVERRIDE_BIN` env var — `runsnake` re-prepends it onto `$PATH`
*after* the module load, so it actually takes priority:

```bash
curl -L -o ~/miniforge3.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash ~/miniforge3.sh -b -p ~/miniforge3
export CONDA_OVERRIDE_BIN="$HOME/miniforge3/bin"
./runsnake 60 -n
```

This doesn't touch the shared/module conda install; `CONDA_OVERRIDE_BIN`
is unset by default and `runsnake` behaves exactly as before if you never
set it.

## Environment / module policy

RepeatMasker, RepeatModeler2, cd-hit, and R all run via **conda or
Singularity only — no `envmodules:` directives anywhere in this
pipeline**. This mirrors the sibling repo `compare_assemblies_satellites`'s
own (deliberate) choice: on this cluster, SGE/DRMAA jobs are submitted with
`-V`, which exports environment *variables* but not shell *functions* — so
`module load` inside a submitted job silently no-ops instead of erroring,
and the job then runs against whatever wrong/absent tool version was
already on `$PATH`. Since nobody has confirmed whether RepeatMasker,
RepeatModeler, or R actually exist as modules on this cluster (`module
avail`), defaulting to the tool stack that's already known to work avoids
that risk entirely. If you confirm real modules exist and want to use them
instead, add `envmodules:` blocks (with the full bootstrap chain below) to
the relevant rules yourself.

If any rule in a *future* version of this pipeline does add
`envmodules:`, the module-load-silently-no-ops issue above requires
loading a **bootstrap chain**, in this order, before any real tool module:

```python
envmodules:
    "modules",
    "modules-init",
    "modules-gs/prod",
    "modules-eichler/prod",
    "sometool/1.2.3",
```

## RepeatModeler2 container

`repeatmodeler.container` in `config.yaml` defaults to
`docker://dfam/tetools:latest`, run via Singularity with `--bind
/net/:/net/` (baked into `runsnake`). This is the only container image used
in the pipeline; `build_db` also runs inside it (not the conda RepeatMasker
env) so the `BuildDatabase`-written database files can't drift to a
different RepeatModeler/RepeatMasker suite version than the one that will
actually read them. Pin this to a specific tag once a real run has checked
which Dfam/RepeatClassifier partitions the image actually bundles — see
`results/{species}/repeatmodeler/{species}.repeatmodeler_provenance.txt`
and `results/summary/provenance.txt`. The tetools image may ship only the
root FamDB partition; if Chordata/Vertebrata content is missing, that's
recorded there, not silently assumed.

**Note on `repeatmodeler`'s `hrs`:** the wiring test is a good place to try
dropping `-l h_rt` for just the `repeatmodeler` job to see if the queue
tolerates an unbounded job. The `runsnake` wrapper submits every job
through one global `--drmaa-args` string, so doing this for one rule only
needs a second, repeatmodeler-specific submission path (not implemented
here) rather than a config toggle. `config.yaml`'s default (`hrs: 96`) is
a normal capped submission; treat any "no h_rt" experiment as a manual,
documented one-off, and record what you find here.

## FamDB / Dfam setup gotchas (from `compare_assemblies_satellites`)

Reused verbatim from that repo's `setup_repeatmasker_famdb` rule (see
`docs/famdb.md` there for the full walkthrough):

1. `FAMDB_DATA_DIR` as a shell/exported env var is **silently ignored** by
   `famdb.py`. The only mechanism that works is a `famdb.conf` file placed
   **next to `famdb.py` itself** — not the working directory.
2. The `-i` flag belongs to **`famdb.py`**, not to `RepeatMasker` — there
   is no supported way to redirect RepeatMasker's own internal FamDB
   lookup from its CLI.
3. Bioconda's RepeatMasker build may bundle its own internal copy of
   `famdb.py`, separate from whatever resolves on `$PATH` — `setup_famdb`
   writes `famdb.conf` next to **every** `famdb.py` found under
   `$CONDA_PREFIX`, not just the one on `$PATH`.
4. `RepeatMasker -species`/`-lib` mode both fail to write a proper `.out`
   header at all if FamDB isn't configured, regardless of whether the run
   actually needs FamDB data (`-lib` mode doesn't) — so `setup_famdb` gates
   every RepeatMasker-invoking rule.

Config keys: `library.famdb_local_dir` (existing local `*.h5` directory) or
`library.famdb_fallback_urls` (download URLs) — confirm the current Dfam
release at https://www.dfam.org/releases/ before filling either in.
`library.dfam_taxon` (default `Vertebrata`) controls the optional
supplementary export (`export_dfam` rule) — expect sparse cyclostome
content; treat it as a supplement, not a foundation, for this lineage.

## `--configfile` target-order caveat

See Quickstart above — this is the same gotcha documented in
`compare_assemblies_satellites`'s README.

## Threading summary

| step | parallelism |
|---|---|
| `BuildDatabase` | none |
| `RepeatModeler -threads` | partial; plateaus (RepeatScout/RECON serial), `-LTRStruct` scales poorly |
| `cd-hit-est -T` | yes |
| `RepeatMasker -pa` | yes; each slot ~4 cores under RMBlast (`repeatmasker.cores_per_pa`) |
| `calcDivergenceFromAlign.pl`, `createRepeatLandscape.pl` | none; run per genome as separate jobs |
| `summarize_rm.py`, plotting | none |
| species-level independence | all per-species rules run concurrently across species |

## Resource tags

Every rule uses the exact pattern:

```python
threads: config["resources"]["<rule>"]["threads"]
resources:
    mem=lambda wildcards, attempt: config["resources"]["<rule>"]["mem"] * attempt,
    hrs=config["resources"]["<rule>"]["hrs"],
    shell_exec="bash"
```

`mem` is **GB per slot** — SGE multiplies by `-pe serial {threads}`
internally, so total memory requested = `mem * threads`. `--restart-times`
must be > 0 (set to 3 in `runsnake`) or the `attempt`-based memory scaling
never engages.

## Queue `h_rt` assumption

`config.yaml`'s resource block assumes a queue that allows at least 96h for
`repeatmodeler` and 48h for `repeatmasker`. **Confirm the target queue's
actual maximum `h_rt` before a real run** — RepeatModeler on a multi-Gb
genome may need more than 96h; see the `hrs` note under "RepeatModeler2
container" above for the "no h_rt" experiment this pipeline hasn't yet
resolved.

## Manual curation hook

`library.curated_override` (path to a FASTA) replaces `shared_library.fa`'s
content entirely with that file (recorded in `provenance.txt`) — use this
after inspecting `results/library/shared_library.fa` and
`results/library/library_membership.tsv` from the `library_only` target,
if some families need manual extension/TSD-checking/TEtrimmer before the
real masking run. `cluster_library` and `library_membership.tsv` still run
unconditionally either way, since the de novo shared-vocabulary comparison
is a result in its own right. TEtrimmer/DeepTE are not implemented in v1.

## Known deviation from the spec's exact file list

`workflow/scripts/` includes one script beyond the spec's original list:
`combine_summaries.py`, which concatenates each `summarize_rm.py` per-arm/
per-species chunk into the final `results/summary/*.tsv` tables and
computes `arm_concordance.tsv`. This keeps `summarize_rm.py` itself focused
on one (arm, species) at a time (matching the spec's description of what
it does) rather than overloading it with cross-run aggregation.

## Known corrections applied vs. the spec's exact rule text (see plan for full rationale)

- `cd-hit-est` is run with `-r 1` (both strands) — the spec's command
  omitted it, which would silently under-merge reverse-complement
  duplicates between two independently-run species (RepeatModeler2 gives
  no guarantee two species' assemblies will report the same family on the
  same strand).
- `build_db` runs inside the same Singularity image as `repeatmodeler`,
  not the conda RepeatMasker env, to avoid a version mismatch between the
  database writer and reader.
- `library_membership.tsv` additionally reports per-cluster
  `label_agreement`/`distinct_labels`, since cd-hit-est keeps one arbitrary
  representative per cluster and silently discards the rest — if two
  species' independent RepeatClassifier calls disagreed on a merged
  family's `Class/Family`, that's now visible instead of silently biasing
  `class_composition.tsv` toward whichever label cd-hit happened to keep.

## Provenance

`results/summary/provenance.txt` records RepeatMasker/RepeatModeler
versions, the FamDB release info, cd-hit version, which Dfam
partition(s) each species' RepeatModeler container could see, any
`curated_override` in effect, and a full `config.yaml` snapshot.
