# Addendum 1 (revised): satellite screen and LTR discovery outside RepeatModeler

This replaces the first draft of Addendum 1 to `repeat_compare_spec.md`. It
records the design as implemented and the decisions made while reviewing
the draft. The README has the full "how"; this document keeps the "what and
why" in one place.

## Why
On *E. stoutii*, RepeatModeler's `-LTRStruct` step (one whole-genome,
single-threaded `gt ltrharvest`) ran for days before it finished. An
independent RepeatMasker pass with KM's satellite library found ~17.9% of
the 2.47 Gb assembly is satellite, spread evenly across every
chromosome-scale scaffold (16.6–18.9% each). That's consistent with
LTRharvest's known slowdown on dense tandem repeats.

Separately, the satellite BED had first been computed against the
scaffolded assembly while RepeatModeler had run on the pre-scaffold (`ptg*`)
assembly. The pipeline must not let that pass silently.

## Decisions (as implemented)

| Topic | Draft said | Implemented |
|---|---|---|
| Canonical assembly | either, but consistent | **scaffolded**, for both species; RepeatModeler is rerun on it (no reuse of earlier runs; `outdir: results_v2`) |
| `-LTRStruct` | config-gated `internal` / `standalone` per species, chosen from satellite % at DAG-build time | **always standalone** for every species; the `repeatmodeler.ltrstruct` key is removed. There's no ordering problem and no checkpoint |
| Masking before LTR discovery | soft-mask, hard-mask as fallback, gated by `satellite_pct_threshold` | **always hard-mask (N)** when a satellite library is set; no threshold. `gt suffixerator` ignores case, and N is its only separator |
| Chunking | whole contigs split into `n_chunks` by bp; `gt suffixerator`/`ltrharvest` per chunk | two levels: `n_groups` (default 8) bp-balanced groups of whole scaffolds as SGE jobs; inside each group, **LTR_HARVEST_parallel** / **LTR_FINDER_parallel** use 5 Mb windows, 100 kb overlap and per-window timeouts. Satellite content is uniform, so a chromosome-sized chunk looks like the whole genome, and only fixed-size windows bound the work per job |
| Which RepeatModeler steps are reused | concatenate `consensi.fa.classified` + `*.mod.LTRlib.fa` | only ltrharvest is replaced. RepeatModeler's own LTRPipeline code runs LTR_retriever, MAFFT, NINJA and Refiner (`LTRPipeline_from_scn`). RepeatModeler's cd-hit merge is ported (`merge_families.py`), and RepeatClassifier runs on the merged set, the same order RepeatModeler uses |
| LTR_retriever input | merged GFF3 | LTRharvest-format `.scn` (what `-inharvest` takes); candidates remapped to the whole-genome 0-based seq-nr LTRPipeline uses |
| ltrharvest parameters | unspecified | RepeatModeler's own (the `gt` defaults); LTR_FINDER added for sensitivity (the one deliberate addition beyond RepeatModeler) |
| Window timeouts | per-chunk `hrs` plus restarts | per-window `-time`. A timed-out window is re-split into 50 kb pieces, and pieces that still time out are **skipped and logged**; the union of skipped bp is reported per species |
| Genome identity guard | md5 of the same FASTA path across three rules | that check can never fail inside one DAG, so it's replaced by a **fingerprint** stored beside the `RM_*` directory and checked before any `-recoverDir` |
| RepeatModeler restart | `-recoverDir` | handles RepeatModeler's exit-0 no-op recover when all rounds are done; refuses a directory built on a different genome |
| Satellite screen | `RepeatMasker -lib … -xsmall -gff -a` | adds `-nolow` (matches the 17.9% estimate; otherwise simple repeats count as satellite); no `-s`; reuses the masking arms' chunks |
| Gap / coverage math | `seqtk cutN` + `bedtools merge` via envmodules | stdlib Python. The pipeline has no envmodules, and `runsnake` doesn't pass `--use-envmodules` |
| Per-contig flag | 40% | **50%**, meaning "mostly satellite"; informational only. The genome-wide % is on every row |
| Satellite library placement | shared library only | **both arms**, like Dfam, so `arm_concordance.tsv` isn't confounded by satellites being in one arm only |
| Duplicate satellite families | accepted limitation | **reclassify, don't cluster**: an `rnd-*#Unknown` family covered by a satellite monomer ×3 target is relabeled `#Satellite`, with the matching harmonized name in a TSV |
| `library_membership.tsv` `source` | tag satellite clusters | satellite entries are never clustered, so it gets a `sources` column (`rnd` / `ltr` / `sat_relabeled`) plus appended rows for Dfam and satellite entries |
| Satellite report | one number | `satellite_composition.tsv` gives an upper bound (satellite-only screen) and a lower bound (shared arm's Satellite class), with `tissue` |
| Species codes | a `code_map` in config | codes read from stage 02b's own `taxonomy_cache/{code}.json` (the collision-resolved code), matched on `species_name`. It's an error if a species is missing, unless `allow_unharmonized` is set. Codes become the library prefix. No NCBI queries |
| 02b outputs | read in place | **snapshotted** into `{outdir}/satellite_harmonization/` with md5s; a changed source re-triggers downstream in one run |
| Per-species override libraries | allowed | allowed, but they go into that species' own arm only, never the shared library |
| RepeatModeler sampling depth | not discussed | defaults for now; `discovery_round_saturation.tsv` shows whether the last round still matters (then use `-numAddlRounds`) |
| Splitting contigs at Ns before RepeatModeler | question raised | **no**: RepeatModeler samples ~40 kb blocks uniformly by non-N bp, so scaffold length doesn't bias sampling |

## Corrections to the draft's factual claims
- The base manifest parser didn't pad to 6 columns; it rejected anything
  other than 5. It now accepts 5 or 6.
- `LTR_retriever` and `gt` exist in `dfam/tetools` but aren't on `PATH`
  (`/opt/LTR_retriever`, `/opt/genometools/bin`). The rules resolve them
  through RepeatModeler's `RepModelConfig.pm`.
- A run killed during `-LTRStruct` has no `consensi.fa.classified`: in
  RepeatModeler, RepeatClassifier runs after the LTR pipeline.
- `MYX2` isn't necessarily genus *Myxine*. Promoted codes carry
  collision suffixes, and a motif shared by *Eptatretus* and *Myxine*
  (two genera) would be promoted to family Myxinidae. Check
  `taxon_codes.tsv` and `harmonized_summary.tsv`; never infer meaning from
  the prefix.
- The upstream LTR_HARVEST_parallel / LTR_FINDER_parallel salvage mode
  re-ran a timed-out window with **no** timeout (its `-threads 1` branch
  ignores `-size` and `-time`). The vendored copies are patched; see
  `workflow/vendor/README.md`.

## Still to confirm on the cluster
- `singularity exec docker://dfam/tetools:2.00 ls /opt/LTR_retriever/LTR_retriever /opt/genometools/bin/gt`
  and that `workflow/scripts/rm_config_path.sh GENOMETOOLS_DIR` resolves
  inside the image.
- The conda `ltr_finder` env's perl has ithreads (the rule checks this and
  fails with a clear message if not).
- A wiring test (`genome_prep.test_subsample_bp` ~50 Mb per species, big
  enough to contain LTRs) through `library_only`. Check that:
  - `rawLTR.scn` passes validation;
  - `ltrs.fa` is non-empty;
  - the `merge_families` log shows RepeatModeler-style counts;
  - `{species}-families.fa` carries RepeatClassifier labels.
- The per-group `hrs` and `threads`, after the first real run's wall-clock
  spread.
