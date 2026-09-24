# repeat_compare

RepeatModeler2 + RepeatMasker cross-species repeat comparison, built for
*Eptatretus stoutii* (de novo) vs *Myxine limosa* (published), but
species-agnostic and manifest-driven.

## Design

Both genomes are masked with **two arms**:

| arm | library | role |
|---|---|---|
| `shared` | non-redundant union of all species' de novo families (+ optional Dfam export + the shared satellite library) | **primary comparison** |
| `own` | that species' de novo families only (+ same optional Dfam export + that species' satellite library) | sanity check / concordance |

Masking every species with only its own library biases the comparison
(each genome is best-annotated for its own families), so `shared` is the
one to trust for cross-species numbers; `own` exists to sanity-check it via
`{outdir}/summary/arm_concordance.tsv`. (`outdir` is `results_v2` by
default: the restructured pipeline below writes to a fresh directory so
runs made with the previous `-LTRStruct` version in `results/` stay
untouched. Paths below are written as `{outdir}/...`.)

Hagfish undergo programmed germline-to-soma genome rearrangement — a
germline assembly and a somatic assembly are different genomes. This
pipeline does **not** act on the manifest's `tissue` column; it only
carries it into every summary table and warns (never fails) if species
disagree or are `unknown`. Assembly-quality covariates (contig count,
total/N/non-N length, N50) are reported alongside every repeat result for
the same reason — a fragmented or collapsed assembly undercounts repeats.

## Pipeline flow

Per species unless noted:

```
prep_genome -> genome_fingerprint, assembly_stats
build_db -> repeatmodeler (RECON/RepeatScout rounds only, no -LTRStruct)
                     |  rounds.consensi.fa / rounds.families.stk (unclassified)
satellite screen (if a satellite library):  split_genome chunks -> satellite_chunk
    (RepeatMasker -nolow -lib <satellite lib>) -> gather_satellite -> satellite_coverage
LTR side pipeline:  ltr_mask (hard-mask satellite hits) -> ltr_group_genome
    -> ltr_harvest_group + ltr_finder_group (per group) -> ltr_gather -> ltr_pipeline
merge_families (rounds + LTR, RepeatModeler's own cd-hit merge)
    -> classify_families (RepeatClassifier) -> satellite_relabel (if a satellite library)
prefix_library -> cluster_library (all species) -> shared / own libraries (+ Dfam, + satellite)
-> repeatmasker (both arms) -> divergence -> summarize -> combine_summaries -> plot
```

### Why LTR discovery runs outside RepeatModeler

RepeatModeler 2.0.9's `-LTRStruct` runs `LTRPipeline`:
1. **One** whole-genome `gt suffixerator` + `gt ltrharvest` with default
   parameters, single-threaded. On the satellite-rich *E. stoutii*
   assembly (~18% satellite, spread evenly over every chromosome-scale
   scaffold) this step ran for days with no output.
2. LTR_retriever (sequences renamed `seqN`, `-noanno`).
3. MAFFT, NINJA and `Refiner` to build `ltr-1_family-N` consensi and seed
   alignments.
4. RepeatModeler then merges those with the round families
   (`cd-hit-est -aS 0.8 -c 0.8 -g 1 -G 0 -A 80`; in a mixed cluster the
   LTR family wins) and runs RepeatClassifier on the merged set.

This pipeline replaces **only step 1**. `ltr_pipeline` runs
`workflow/vendor/RepeatModeler/LTRPipeline_from_scn` (RepeatModeler's own
`LTRPipeline`, patched to read a precomputed `.scn`) for steps 2–3.
`merge_families.py` ports step 4's merge, and `classify_families` runs
RepeatClassifier. Every species goes through the same path, so results are
comparable across species. Relative to a stock `-LTRStruct` run, two things
differ, and methods should say so:
- LTR candidates come from the **satellite-hard-masked** genome;
- **LTR_FINDER** candidates are added to LTRharvest's
  (`ltr_discovery.use_ltr_finder`).

`ltr_discovery.ltrharvest_args: ""` keeps RepeatModeler's own ltrharvest
parameters (the `gt` defaults). `merge_families.py` also fixes two
RepeatModeler edge cases, which its docstring lists: the last cd-hit
cluster was never evaluated, and LTR families were dropped when no family
was redundant.

Parallelism has two levels:
- **SGE jobs:** `ltr_group_genome` splits *whole* scaffolds into
  `ltr_discovery.n_groups` bp-balanced groups (default 8). Scaffolds are
  never split across groups, so merging groups is a plain concatenation.
- **Threads within a job:** the vendored LTR_HARVEST_parallel and
  LTR_FINDER_parallel cut each group into 5 Mb windows with 100 kb
  overlap, and kill any window after `window_timeout_s`. With `try1: 1`, a
  killed window is re-run as 50 kb pieces with their own timeouts. Pieces
  that still time out are skipped and logged.

`ltr_gather` (`normalize_scn.py`) rewrites every candidate to the 0-based
whole-genome sequence index LTRPipeline expects, validates every
coordinate against the sequence length, and reports the union of skipped
sequence. The result goes to `{outdir}/summary/ltr_discovery.tsv` and
`provenance.txt`. See `workflow/vendor/README.md` for what was patched in
the vendored tools. That includes an upstream bug: salvage mode re-ran a
timed-out window with no timeout at all.

**If LTR candidate jobs are slow or failing, check these first:**
- `{outdir}/{species}/ltr/groups/manifest.tsv` shows which contigs are in
  which group;
- `{outdir}/{species}/ltr/skipped_windows.tsv` shows which regions timed
  out;
- the per-group `*.timeouts.tsv` files list every salvaged or skipped piece.

Tool paths: in `dfam/tetools`, `gt`, `LTR_retriever` and `cd-hit` live
under `/opt` but are **not on `PATH`**, so `LTR_retriever` on its own says
"command not found". The rules resolve them the way RepeatModeler does,
through its `RepModelConfig.pm` (`workflow/scripts/rm_config_path.sh`).
`ltr_finder` comes from bioconda (`workflow/envs/ltr_finder.yaml`).

### Satellite screen, and the satellite library in the masking libraries

Species with a satellite library (see below) get a satellite-only
RepeatMasker screen: `-nolow -lib <library>`, no `-s`, the settings behind
the original 17.9% estimate. It reuses the `split_genome` chunks. Outputs:
- `{outdir}/{species}/satellite_screen/satellite_genomewide.tsv` and
  `satellite_per_contig.tsv`: satellite bp against total and non-gap bp.
  N runs are subtracted, so no seqtk or bedtools module is needed. The
  per-contig table sets `flagged` at ≥ `satellite.per_contig_flag_pct`
  (50, meaning "mostly satellite", informational only) and repeats the
  genome-wide % on every row.
- `satellite.bed`: merged hits, hard-masked by `ltr_mask` for LTR
  discovery only.
- `{outdir}/summary/satellite_composition.tsv`, with both bounds:
  - `*_screen` (upper bound): the satellite library alone, no competition;
  - `*_shared_arm` (lower bound): the shared arm's `Satellite` class, where
    the same entries compete with de novo and Dfam families. Some satellite
    loci go to an `rnd-*#Unknown` duplicate there.
- `{outdir}/summary/satellite_by_origin.tsv` (with `harmonization_dir`):
  satellite bp split into the species' own motifs, motifs promoted to a
  shared ancestor (with `lca_rank`/`lca_taxon`), and other species'
  motifs.

The satellite library is appended to **both** masking arms, unclustered
and unprefixed, the same way Dfam is. That keeps `arm_concordance.tsv` a
comparison of clustering effects rather than of which arm has satellites.
`append_libraries.py` checks every header is `name#Class/Family`,
deduplicates libraries by content, skips names a `curated_override`
already contains, and fails on name collisions.

The satellite library is **not** cd-hit-clustered against de novo
families, because cd-hit has no notion of repeat phase or rotation.
Instead, `satellite_relabel` screens the classified families against each
satellite monomer concatenated ×3. An `rnd-*#Unknown` family covered at
least `satellite.relabel_min_cov` is relabeled `#Satellite`, and the
matching harmonized name goes into `families/satellite_relabel.tsv`.
Harmonized entries are never renamed. Some fragmentation can remain at
family level, where one satellite appears under both a harmonized name and
an `rnd-*` name. Class-level `Satellite` bp is the comparable number.

### Satellite library QC (does each motif behave like a satellite here?)

`satellite_library_qc` checks every library motif against the operational
definition used by `compare_assemblies_satellites`. It measures each
criterion on this genome's satellite screen. It's **report-only**: nothing
is filtered, and masking uses the full library.

| criterion | measured as | default (`satellite.qc`) |
|---|---|---|
| length | monomer length | 75–2000 bp (stage 02 ran with `min_period_length: 75`; 2000 is TRF's maximum period) |
| copies | genome-wide: merged hit bp / monomer length | ≥ 100 |
| tandem | fraction of the motif's bp in arrays (same motif, contig and strand, hits chained across gaps ≤ max(50 bp, 0.2 × monomer)) spanning ≥ 3 monomers | ≥ 0.5 |
| not a simple repeat | fraction of the monomer covered by an exact period-1..10 self-repeat (≥ 3 copies of the unit) | ≤ 0.5 |

Stage 02's copy cutoffs (`candidate_scan.min_copy_number: 100`,
`min_single_block_copy_number: 300`, `copy_number_filter.min_total_copy_number: 500`)
are TRF copy numbers **within arrays**. They don't establish genome-wide
copy number or tandem organisation for a harmonized motif screened against
a genome, which is why this recomputes both.
- A motif whose bp is mostly isolated hits behaves like a dispersed
  repeat, often a TE fragment. It still inflates satellite %, and it gets
  hard-masked before LTR discovery.
- A monomer that is mostly a short-period repeat is effectively
  `Simple_repeat`. The screen's `-nolow` lets it count as satellite.

**Is 100 copies enough?** As an inclusion floor, yes. As a definition of
"a real satellite", no:
- 100 copies of a 75–2000 bp monomer is only 7.5–200 kb, a few thousandths
  of a percent of a 2.5 Gb assembly. Major satellites run to 10⁴–10⁶
  copies in kb–Mb arrays.
- Any motif stage 02 discovered in a genome clears 100 copies in that same
  genome almost automatically. The floor mainly catches harmonized motifs
  borrowed from another species that barely occur here.
- A fixed copy count favours long monomers.

So pass/fail stays at `min_copies` and each motif also gets a `copy_tier`:
- `major`: ≥ `major_min_copies` (1000) copies, **or** ≥ `major_min_bp`
  (100 kb) of sequence;
- `minor`: passes `min_copies` but isn't major;
- `below_floor`: under `min_copies`.

A long monomer can be `major` by bp while failing `pass_copies` (80
copies of 1.8 kb is 144 kb). That's the long-monomer effect made visible,
not a bug.

`plots/satellite_copy_distribution.png` plots copies against monomer
length for every motif, one panel per species, with each cutoff drawn in.
If the copies show a clear gap, move `min_copies` there. That's the same
idea as stage 02's `copy_number_diagnostic.png`.

**Outputs:**
- Early: `{outdir}/{species}/satellite_screen/satellite_library_qc.tsv`
  (one row per motif, a `pass_*` flag per criterion, `pass_all`,
  `copy_tier` and `largest_array_bp`).
- Final: `{outdir}/summary/satellite_library_qc.tsv`, which adds
  `te_like_family_hits`: de novo families that RepeatClassifier calls a TE
  and whose best satellite match is this motif, a hint the motif is
  TE-derived. This needs RepeatModeler, so it only appears in the final
  summary.
- `satellite_composition.tsv` gains two columns:
  `satellite_pct_nongap_screen_passing` (only motifs that pass every
  check) and `satellite_pct_nongap_screen_major` (only major-tier motifs).
  Together they show how much the headline number depends on borderline
  and minor motifs.

**Check it before the full run:** the `satellite_qc` target runs only
prep → satellite screen → QC. That takes hours, compared with days for
RepeatModeler:

```bash
./runsnake 40 --configfile config.yaml satellite_qc
```

### Species codes and stage 02b harmonization

Set `satellite.harmonization_dir` to
`compare_assemblies_satellites/results/02b_taxonomy_harmonization`. Then:
- `snapshot_harmonization` copies its harmonized library, summary, code
  tables, taxonomy tree and `taxonomy_cache/` into
  `{outdir}/satellite_harmonization/`, with an md5 manifest. Every later
  rule reads only that copy. Changing a source file re-triggers the
  snapshot and everything downstream in a **single** run.
- Each species' code (EST, MLI, ...) comes from 02b's
  `taxonomy_cache/{code}.json`, by matching the manifest `species_name`
  against the binomial 02b was given. The file name is the
  collision-resolved code 02b actually used. The code becomes the library
  family prefix (`EST_rnd-1_family-3`, matching `EST_SAT...`) and the
  `species` column of every table. Paths keep `species_id`.
- A species that isn't found is an error listing the available binomials.
  Set `satellite.allow_unharmonized: true` to fall back to `species_id`
  instead.
- Without `harmonization_dir`, codes are the manifest `species_id`s and
  satellite libraries come from manifest column 6 only.

Never read meaning into a code prefix. Promoted codes are the first three
letters of the ancestor's name plus a collision suffix, so `MYX` / `MYX2`
can be genus *Myxine* or family Myxinidae depending on collision order.
`satellite_by_origin.tsv` takes `lca_rank`/`lca_taxon` from
`harmonized_summary.tsv` instead. No NCBI lookups happen here: 02b's cached
lineages produced the names and are the source of truth.

A manifest column-6 library is a per-species **override**. It goes only
into that species' own arm and satellite screen, never into the shared
library, so one species' un-harmonized motifs never mask every species.

### Genome identity guard and RepeatModeler restarts

`genome_fingerprint` records the prepped FASTA's md5 and every sequence
name and length. Before starting a run, `repeatmodeler` stores that
fingerprint next to the `RM_*` directory (`rm_run.fingerprint.tsv`). On any
later attempt it refuses to touch an existing `RM_*` directory unless the
current genome matches. Without this, a changed input (pre-scaffold vs
scaffolded assembly, or a changed `genome_prep.test_subsample_bp`) would be
silently `-recoverDir`-ed on top of a run built from a different genome.

RepeatModeler's `-recoverDir` exits 0 without doing anything when every
round already finished ("appears to contain a successful run"). The rule
accepts that case, since only the rounds are needed, and otherwise
requires a real completion.

**Rule of thumb:** the satellite screen, RepeatModeler and the LTR side
pipeline must all run on the identical genome FASTA for a species, every
time. Scaffolded or pre-scaffold doesn't matter, as long as it's the same
file. The pipeline guarantees this within a run. The guard catches it
across runs.

### RepeatModeler sampling depth

RepeatModeler doesn't sample whole sequences. It cuts every sequence into
~40 kb blocks, shuffles all blocks together, and draws blocks without
replacement until each round's **non-N** bp target is reached (RepeatScout
40 Mb, then RECON 10 → 30 → 90 → 270 Mb, about 15–20% of a 2.5 Gb
assembly). Every region is therefore sampled in proportion to its non-N
bp, whatever its scaffold's length.
- Unplaced contigs aren't under-sampled. One shorter than 40 kb is a
  single block, so it's slightly *over*-represented per bp.
- Splitting contigs at Ns first would only dilute them.

`{outdir}/summary/discovery_round_saturation.tsv` shows own-arm masked bp
by the round that discovered each family (`rnd-1` … `rnd-N`, `ltr`,
`satellite`, `other`). If families from the final round still mask a
meaningful share (for example >1% of the genome), sampling hasn't
saturated. In that case set `repeatmodeler.extra_args: "-numAddlRounds 1"`
(or 2), the same value for every species, and rerun.
- Each extra round costs about as much as the most expensive round, and
  later rounds mostly add low-copy `Unknown` families.
- Prefer extra rounds over a larger `-genomeSampleSizeMax`: RECON's
  all-vs-all cost grows superlinearly with sample size.

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

`repeatmodeler.container` in `config.yaml` is pinned to
`docker://dfam/tetools:2.00`, the tag whose digest matched `latest` for the
first runs (RepeatModeler 2.0.9, LTR_retriever 2.9.0, genometools 1.6.4).
It runs via Singularity with `--bind /net/:/net/`, which is baked into
`runsnake`. This is the only container image in the pipeline.
`build_db`, `repeatmodeler`, the LTR rules (`ltr_harvest_group`,
`ltr_pipeline`), `merge_families` and `classify_families` all run inside
it, so no step drifts to a different RepeatModeler/RepeatMasker suite
version. To check which Dfam/RepeatClassifier partitions the image
actually bundles, see
`{outdir}/{species}/repeatmodeler/{species}.repeatmodeler_provenance.txt`
and `{outdir}/summary/provenance.txt`. The tetools image may ship only the
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

## RepeatMasker scatter/gather

Each species' genome is masked via a 3-rule scatter/gather
(`split_genome` → `repeatmasker_chunk` × `repeatmasker.scatter_count` →
`gather_repeatmasker`), not as one unchunked whole-genome job — adapted
from `compare_assemblies_satellites`'s stage 03 `-lib`-mode RepeatMasker
scatter/gather, which solves the identical scheduling problem: many
small, independently-restartable SGE jobs schedule and recover from
failure much better than one large reservation for the whole genome.
`workflow/scripts/split_fasta.py` and `workflow/scripts/gather_rm_out.sh`
are copied verbatim from that repo's `common/scripts/`. One thing
deliberately **not** carried over: that stage's `-nolow` flag (it
suppresses RepeatMasker's built-in low-complexity/simple-repeat screen) —
we need that screen to run, since `Simple_repeat` and `Low_complexity`
are both required canonical output classes here.

`gather_repeatmasker` merges three file types per (arm, species):
`.out` via `gather_rm_out.sh` (keeps one header, strips the zero-hit
sentinel line per chunk); `.tbl` by summing each chunk's "total
length"/"bases masked" numbers into a minimal synthetic file (chunks are
disjoint contig sets, so these are additive, and `summarize_rm.py`'s
`.tbl` cross-check only ever greps for "bases masked" anyway); `.align`
by straight concatenation. The `.align` concatenation has no proven
precedent in the sibling repo (it never needed to merge that file type) —
verify it produces valid `calcDivergenceFromAlign.pl` input during the
wiring test before trusting it for the full run.

`repeatmasker.scatter_count` (default 10) and the `repeatmasker` resource
block are **per chunk now**, not per whole genome — untuned placeholders,
adjust both from real per-chunk runtimes observed in the wiring test.

## Threading summary

| step | parallelism |
|---|---|
| `BuildDatabase` | none |
| `RepeatModeler -threads` | partial; plateaus (RepeatScout/RECON serial); rounds only (no `-LTRStruct`) |
| satellite screen (`RepeatMasker -pa`) | yes, per chunk, like the masking arms |
| `ltr_harvest_group` / `ltr_finder_group` | yes: `n_groups` SGE jobs per species × `threads` windows each (threads must be ≥ 2) |
| `ltr_pipeline` (LTR_retriever, MAFFT, NINJA) | yes, `-threads`; once per species (needs the whole genome's candidates) |
| `merge_families` (cd-hit-est), `classify_families` (RepeatClassifier) | yes, `-T` / `-threads` |
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
after inspecting `{outdir}/library/shared_library.fa` and
`{outdir}/library/library_membership.tsv` from the `library_only` target,
if some families need manual extension/TSD-checking/TEtrimmer before the
real masking run. `cluster_library` and `library_membership.tsv` still run
unconditionally either way, since the de novo shared-vocabulary comparison
is a result in its own right. TEtrimmer/DeepTE are not implemented in v1.

## Known deviation from the spec's exact file list

`workflow/scripts/` includes one script beyond the spec's original list:
`combine_summaries.py`, which concatenates each `summarize_rm.py` per-arm/
per-species chunk into the final `{outdir}/summary/*.tsv` tables and
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

`{outdir}/summary/provenance.txt` records RepeatMasker/RepeatModeler
versions, the FamDB release info, cd-hit version, which Dfam
partition(s) each species' RepeatModeler container could see, the pinned
container tag, each species' genome fingerprint, the LTR tool versions
(genometools, LTR_retriever, ltr_finder; vendored script commits are in
`workflow/vendor/README.md`), LTR candidate counts and window-timeout
skipped bp, the species code map, which satellite library went into which
arm (plus the stage-02b snapshot md5s), any `curated_override` in effect,
and a full `config.yaml` snapshot.
