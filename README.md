# repeat_compare

RepeatModeler2 + RepeatMasker cross-species repeat comparison, built for
*Eptatretus stoutii* (de novo) vs *Myxine limosa* (published), but
species-agnostic and manifest-driven.

## Design

Both genomes are masked with **two arms**:

| arm | library | role |
|---|---|---|
| `shared` | non-redundant union of all species' de novo families (+ optional Dfam export) | **primary comparison** |
| `own` | that species' de novo families only (+ same optional Dfam export) | sanity check / concordance |

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
LTR side pipeline:  ltr_group_genome
    -> ltr_harvest_group + ltr_finder_group (per group) -> ltr_gather -> ltr_pipeline
merge_families (rounds + LTR, RepeatModeler's own cd-hit merge)
    -> classify_families (RepeatClassifier)
prefix_library -> cluster_library (all species) -> shared / own libraries (+ Dfam)
-> repeatmasker (both arms) -> divergence -> summarize -> combine_summaries -> plot
```

### Why LTR discovery runs outside RepeatModeler

RepeatModeler 2.0.9's `-LTRStruct` runs `LTRPipeline`:
1. **One** whole-genome `gt suffixerator` + `gt ltrharvest` with default
   parameters, single-threaded. On the highly repetitive *E. stoutii*
   assembly this step ran for days with no output.
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
- LTRharvest runs on **overlapping 5 Mb windows with per-window timeouts**
  instead of one whole-genome pass. Windows that still time out are
  skipped, and their bp are reported in `{outdir}/summary/ltr_discovery.tsv`.
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

**Rule of thumb:** RepeatModeler and the LTR side pipeline must both run
on the identical genome FASTA for a species, every
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
`other`, where `other` means Dfam and RepeatMasker's own simple/low-complexity
calls). If families from the final round still mask a
meaningful share (for example >1% of the genome), sampling hasn't
saturated. In that case set `repeatmodeler.extra_args: "-numAddlRounds 1"`
(or 2), the same value for every species, and rerun.
- Each extra round costs about as much as the most expensive round, and
  later rounds mostly add low-copy `Unknown` families.
- Prefer extra rounds over a larger `-genomeSampleSizeMax`: RECON's
  all-vs-all cost grows superlinearly with sample size.

### Discovery summary

`{outdir}/summary/discovery_summary.tsv` follows the family counts from
discovery to the shared library, with one row per species plus `ALL`.
It's built from library-stage files only, so `library_only` produces it too.

| Column | Meaning |
|---|---|
| `rounds_families`, `ltr_families` | families from RepeatModeler's rounds / the LTR side pipeline |
| `merge_removed`, `merged_families` | round families dropped as redundant with an LTR family, and what's left (`merge_families`) |
| `merged_ltr_families`, `putative_subfamilies` | LTR families in the merged set; round families tagged "putative subfamily of" |
| `classified_families`, `unknown_families` | after RepeatClassifier; those labelled `Unknown` |
| `clusters`, `species_only_clusters`, `shared_clusters` | cross-species cd-hit-est clusters (`cluster_library`) containing this species' families |
| `families_in_species_only_clusters`, `families_in_shared_clusters`, `pct_families_in_shared_clusters` | where this species' families landed |
| `shared_clusters_label_conflict` | shared clusters whose members' `Class/Family` labels disagree (see `library_membership.tsv`) |
| `dfam_entries` | Dfam families appended to the libraries, counted by unique name (`ALL` row only) |

In the `ALL` row, family counts are summed over species and cluster counts
are over the whole clustering. A species' family count must equal its
members in the `.clstr`; otherwise the rule fails, since the inputs would
come from different runs.

`discovery_summary_by_class.tsv` splits families into species-only vs
shared clusters per Class (the part of each family's own label before
`/`), for each species and `ALL`. Sharing is usually very uneven across
classes, so read this alongside the overall percentage.

### Overlap with Dfam

RepeatClassifier compares families with the Dfam in the container, but it
only assigns a class. It doesn't record which Dfam family matched, and the
Dfam export is appended to the libraries without being compared with the
de novo families. `dfam_overlap` (run when `library.include_dfam` is set)
fills that gap. It runs `cd-hit-est-2d` of each species' families against
the Dfam export, with `cluster_library`'s identity threshold, on both
strands:
- `{outdir}/summary/dfam_overlap.tsv`: per species and `ALL`, and per
  Class, how many families match a Dfam family, split into
  `class_agrees` / `class_differs` (family classified; Dfam Class the
  same / different) and `unknown_matched` (family was `Unknown`, so the
  match suggests a class);
- `{outdir}/library/dfam_overlap/dfam_matches.tsv`: each matching family
  with its Dfam family, identity, coverage of the de novo family and
  strand.

A match means the alignment covers at least 80% of the **de novo** family
(`library.cdhit.coverage_short`). A family that is a fragment of a known
element counts. A family that only *contains* a short Dfam entry (a tRNA,
MITE, solo LTR or satellite inside a longer element) doesn't. cd-hit's own
`-aS` measures coverage of the *shorter* sequence, which would accept
those, so the rule adds `-s2 0.8` and the script checks the coverage
itself. Expect few matches in lineages Dfam barely covers, such as
cyclostomes.

## Tandem check (`family_tandem`)

Every family in the shared-arm `.out` of each species is checked against
the operational satellite definition the removed satellite arm used
(`satellite_library_qc.py`, tag `satellite-arm-v1`), so a shared-library
family and a satellite-library motif are judged the same way:

| criterion | rule (config `family_tandem`) |
|---|---|
| length | consensus 75–2000 bp |
| copies | ≥ 100 genome-wide (merged hit bp / consensus length) |
| tandem | ≥ 50% of the family's bp in arrays: same-family hits on one contig and strand, chained when the gap is ≤ max(50 bp, 0.2 × consensus), spanning ≥ 3 consensus lengths |
| not simple | ≤ 50% of the consensus is a period-1..10 exact self-repeat |

`satellite_like` = all four. Fragments of one interrupted TE chain too, but
span about one consensus length, so they never reach 3 copies.

`tandem_family` = copies + tandem + not simple at **any** consensus length.
The 2000 bp cap is TRF's period limit, not biology: on *M. limosa* the
largest satellite, `Mlim_rnd-1_family-3` (2.9 kb consensus, 542 Mb, 99% in
arrays up to 1.85 Mb), fails only that test. With
`family_tandem.carve_unknown: true` (default), `summarize` reports Unknown
families with `tandem_family` as their own class, **Unknown_tandem**
("Unknown (tandem)", light violet), in the composition, landscape and
concordance tables and plots. `family_tandem` runs on both arms for this,
each arm's own table deciding its families. It is a placeholder until the
satellite library's RepeatMasker runs are folded in. `monomer_period` is the
consensus's strongest internal repeat period (k-mer spacing), so a
consensus that is several copies of a shorter unit reports that unit.

- `{outdir}/summary/family_tandem.tsv` (or `summary_shared_only/`): one row
  per family and species, largest first: `owned_bp` (the bp it holds under
  the class table's highest-score rule), tandem fraction, array count and
  largest array, `median_div` (bp-weighted `.out` perc. div., to place a
  family on the landscape), the four pass flags, `satellite_like` and a
  copy tier (major ≥ 1000 copies or ≥ 100 kb).
- `{outdir}/summary/class_tandem.tsv`: per class, how much of the class's bp
  is held by `satellite_like` and by `tandem_family` families, e.g. how much
  of Unknown behaves like satellite.

Apart from the Unknown_tandem class, nothing is relabelled. To check it against the satellite
pipeline, compare `satellite_like` families with the motifs that pass in
that pipeline's own QC on the same assembly.

## Reclassifying Unknown families and verifying classes (`classify`)

Three cheap, independent screens run once on every **de novo** consensus in
the shared library (`{outdir}/classify/`). These are the species-prefixed
RepeatModeler families. Curated Dfam entries keep their labels and are left
out of both the reclassification and the verification;
`classify.screen_dfam: true` includes them. Their share of masked bp is in
`library_source.tsv`. They serve two purposes:
- For **Unknown** families, `reclassify_unknown.py merge` turns the screens
  into `unknown_reclassification.tsv`, which `summarize` applies to the
  `.out`/`.align` labels.
- For **classified** families, `verify` checks RepeatClassifier's labels
  against the same evidence. This is report-only (see "Verifying
  RepeatModeler's classes" below). Family names don't change, so **nothing is
remasked**: turning `classify.enabled` on or off only reruns `summarize` and
everything downstream of it.

| screen | tool | answers | evidence |
|---|---|---|---|
| protein domains | TEsorter, REXdb-metazoa and GyDB, `-dp2` (domain hits only, no similarity pass) | which TE order / superfamily | `domain` |
| structured RNA | Infernal `cmscan --rfam --cut_ga` against Rfam | rRNA / tRNA / snRNA / srpRNA | `rfam` |
| host genes | DIAMOND blastx against TE-cleaned host proteins | is it a gene, not a TE | `host_protein` |

Precedence (first match wins; only `Unknown` families change):

1. **Unknown_tandem** (`family_tandem`, array evidence) always wins.
2. **domain**: a TEsorter call with at least `tesorter_min_domains` domains,
   mapped to RepeatMasker names (`LTR/Gypsy`, `DNA/hAT`, `DNA/CMC-EnSpm`, ...;
   unmapped superfamilies become `<order>/<TEsorter name>`).
3. **rfam**: the hit covers at least `min_rfam_cov` of the consensus.
4. **host_protein**: the best host protein covers at least `min_host_cov`,
   giving `Other/host_gene`. These are excluded from TE classes and counted in Other.

The family stays Unknown, with the reason in `conflict`, if:
- REXdb and GyDB disagree on order, or
- a domain call coincides with a qualifying host or Rfam hit (e.g. a
  domesticated TE gene).

`diamond_best` and `rfam_best` show the best hit even below the coverage cut.
For example, a tRNA-headed SINE shows a partial tRNA hit.

`family_composition.tsv` gains `bp_from_unknown`: the bp of each
Class/Family that came from relabelled Unknown families. Report how much
moved with it, together with the `evidence` column.

**Host proteins.**
- `annotation_proteins` takes the species' gene annotations, e.g. NCBI's
  `complete.proteins.faa` for Esto and the MLI protein FASTA. Not the
  transcripts: DIAMOND needs proteins.
- `swissprot_fasta` is optional.
- Gene annotations still carry TE-derived models, so the proteins are
  filtered before the DIAMOND database is built:
  - all inputs: by `te_protein_keywords` in their descriptions
  - annotation proteins also: by TEsorter's protein mode (any TE domain drops the protein)
- Domesticated TE genes (PGBD, ZBED, CENP-B...) are dropped too, on purpose: a
  family matching one keeps its TE label rather than becoming a "host gene".
- The host log reports how many proteins each filter removed.

**Reference downloads.** Nothing needs installing: `--use-conda` builds the
TEsorter, DIAMOND and Infernal environments. TEsorter's HMM databases ship
with the package. Point the config at the downloads exactly as they come;
`.gz` is fine for every FASTA and for `Rfam.cm`:
- Swiss-Prot: `https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz`
- Rfam: `https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz` (no
  `Rfam.clanin` needed: only each consensus's best hit is used)

`rfam_prep` copies or gunzips `Rfam.cm` into `{outdir}/classify/rfam/` and runs
`cmpress` there, in the pipeline's own Infernal environment. Nothing is
written next to your reference files.

**Notes.**
- A screen whose input is empty (`annotation_proteins`/`swissprot_fasta`, or
  `rfam_cm`) is skipped and treated as "no evidence".
- TEsorter always runs while `classify.enabled` is on.
- The table comes from the shared library. In the own arm it relabels the
  same family names, which covers most families. Own-library families that
  were merged away during clustering keep their original label.
- To toggle from the command line, use e.g.
  `--config "classify={enabled: false}"`. Booleans passed that way arrive as
  strings and are parsed accordingly.

### Verifying RepeatModeler's classes

RepeatClassifier labels families by homology to Dfam and RepeatPeps. The
screens are independent of that, so `verify` compares each classified
family's label with them. It changes nothing; it only reports.
Families are weighted by the bp they hold in each species (`owned_bp` from
the shared-arm `family_tandem.tsv`).

| status | meaning |
|---|---|
| `agree_superfamily` | the domain call matches the label (`LTR/ERV1` vs `LTR/ERV` counts) |
| `agree_order` | same order, different or unspecified superfamily (e.g. `LTR/Gypsy` vs Copia domains) |
| `disagree_order` | the domains point to another order, e.g. a `SINE/Alu` with Gypsy domains |
| `retroposon_vs_line` | a `Retroposon` label with LINE domains: may be an autonomous LINE, or a non-autonomous element carrying LINE fragments |
| `domain_conflict` | REXdb and GyDB disagree on order |
| `host_protein` | no domain, but the consensus matches a TE-free host protein: possible gene family or contamination |
| `rfam` / `agree_rfam` / `rfam_other_rna` | a structured-RNA hit on a non-RNA label / on a matching RNA label / on a different RNA label |
| `domain_on_rna_label` | an RNA label with TE domains |
| `no_evidence` | the screens are silent |

**`no_evidence` is not a failure.** SINEs, MITEs, solo LTRs, satellites
and decayed copies have no protein domains. Read it as "not confirmed",
never as "wrong".

Outputs:
- `{summary}/class_verification.tsv`: per species, class and status, the
  families present, the bp they hold, and that bp as a % of the class.
- `{summary}/class_disagreements.tsv`: families with a disagreeing status
  (`disagree_order`, `retroposon_vs_line`, `domain_conflict`,
  `host_protein`, `rfam`, `rfam_other_rna`, `domain_on_rna_label`), plus
  tandem families with a TE label. Sorted by bp, with each screen's best
  hit. Review this list by hand, starting from the top.

## Library sources (`library_source.tsv`)

`{summary}/library_source.tsv` splits each species' masked bp by where the
family holding it came from:

| source | meaning |
|---|---|
| `own_denovo` | families discovered in this species (`<species_id><sep>...`) |
| `denovo:<species>` | families discovered in the other species |
| `dfam` | the Dfam export (`library.include_dfam` / `dfam_taxon`) |
| `rm_builtin` | RepeatMasker's own simple-repeat / low-complexity screen |
| `total` | all of the above |

Columns:
- `n_families`: families holding any bp
- `n_hits`: `.out` hit lines (copies and fragments), the N for a methods table
- `owned_bp`
- `pct_of_masked`
- `pct_non_n`: bp over the non-N assembly length
- `bp_weighted_median_div`: how old each source's matches are

bp follow the same rule as `class_composition.tsv`: every base counts once,
for its best hit, so the `total` row matches the composition total.

Two caveats:
- Each base counts toward the family that *won* it at masking time. The Dfam
  export is appended to the library after the de novo families are
  clustered, so a Dfam entry and a near-identical de novo family both stay
  in the library and compete base by base. Where the de novo consensus fits
  better, the bp count as de novo even if Dfam has the same element.
  `dfam_overlap.tsv` lists de novo families that match known Dfam families.
- A family shared by both species is one cd-hit cluster, kept under its
  longest member's name, so `denovo:<other species>` includes shared
  families whose representative came from the other species
  (`library_membership.tsv` has the clusters).
- On the first *E. stoutii* / *M. limosa* run, Dfam held under 1% of masked bp
  in both species (old, partial matches, ~24% divergence), while about 11% of
  each genome was masked by the *other* species' de novo families.

## Satellite analysis (removed)

A satellite arm existed briefly. It was a satellite-only RepeatMasker
screen with the harmonized library from `compare_assemblies_satellites`
stage 02b, per-motif QC, satellite masking before LTR discovery, and the
satellite library appended to both masking arms. It was removed because
its first QC run on *E. stoutii* showed most of the library's
satellite-screen bp weren't tandem:
- 181 high-copy motifs, carrying 86% of screen bp, were low-divergence,
  partial, mixed-strand fragments dispersed genome-wide;
- only about 72 Mb, roughly 3% of the genome, sat in tandem arrays,
  against 17.8% from all hits.

The satellite caller will be tightened first. The last commit that
includes the satellite arm is `502accf`, tagged `satellite-arm-v1`
where the tag has been pushed. To bring it back:

```bash
git checkout 502accf -- Snakefile config.yaml workflow/scripts   # or cherry-pick specific files
```

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

**Mask with the shared library only (`mask_shared_only`):** to run just
RepeatMasker with an already-built shared library on every species, with no
discovery, LTR, clustering, own-arm or summary rules, point
`mask_shared_library` at the library and ask for the `mask_shared_only` target:

```bash
./runsnake 40 --configfile config.yaml \
    --config mask_shared_library=results/library/shared_library.fa \
    -- mask_shared_only
```

Everything after `--` is a target, so flags like `-n` must go *before* the
`--` (`... shared_library.fa -n -- mask_shared_only`).
Only `prep_genome`, `split_genome`, `setup_famdb` (quick ones) and the
shared-arm `repeatmasker_chunk`/`gather_repeatmasker` jobs are scheduled.
Output goes to `<outdir>/shared/<species>/repeatmasker/`. Leave
`mask_shared_library` unset for a normal run.

**Report on the shared arm only (`report_shared_only`):** after
`mask_shared_only`, run the reporting that needs only the shared arm, with
the **same** `mask_shared_library` value:

```bash
./runsnake 10 --configfile config.yaml \
    --config mask_shared_library=results/library/shared_library.fa \
    --rerun-triggers mtime -n -- report_shared_only   # check: no repeatmasker_chunk jobs listed
./runsnake 10 --configfile config.yaml \
    --config mask_shared_library=results/library/shared_library.fa \
    --rerun-triggers mtime -- report_shared_only
```

Keep `--rerun-triggers mtime`. The per-chunk RepeatMasker outputs are temp
files, deleted after `gather_repeatmasker`. If Snakemake's recorded metadata
makes it want to redo the gather (for example "Set of input files has
changed since last execution"), it has to re-mask every chunk to rebuild
them. With `mtime` it only reruns the gather if the genome or library is
actually newer than the masked output.

`mtime` also ignores script changes, so after a pipeline update that changes
`summarize_rm.py` (for example the overlap-resolved divergence landscape),
add `--forcerun summarize`. That reruns the per-species summaries and
everything downstream, but not the masking.

The target refuses to run without `mask_shared_library`. Without it the
shared arm would mask with `<outdir>/library/shared_library.fa`, which
schedules the whole discovery chain and re-runs RepeatMasker.

It schedules `assembly_stats`, `divergence` and `summarize` per species,
then `combine_summaries_shared_only` and `plot_shared_only`:
- `<outdir>/summary_shared_only/`: `class_composition.tsv`,
  `family_composition.tsv`, `divergence_landscape.tsv`,
  `assembly_covariates.tsv` (same formats as `summary/`);
- `<outdir>/plots_shared_only/`: `class_composition_shared.png`,
  `divergence_landscape.png`;
- `<outdir>/shared/<species>/divergence/<species>.landscape.html`:
  RepeatMasker's own landscape page.

These go to separate directories so a later full run's `summary/` and
`plots/` are never mixed with them. Not produced: `arm_concordance` (needs
the own arm), `discovery_round_saturation` (own arm + RepeatModeler
rounds), `ltr_discovery`, `discovery_summary`, `dfam_overlap` and
`provenance.txt` (discovery chain).

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

**Jobs inherit the launching shell's `PATH`.** `runsnake` submits with
`-V`, so every job gets the `PATH` of the shell that started Snakemake,
including whether you were in `(base)` or had run `conda deactivate`.
From a `(base)` shell, Miniforge's base `python3` can come before a rule
env's, so `famdb.py` (`#!/usr/bin/env python3`) fails with
`No module named 'h5py'` even though the rule's env has it.

The rules that use `workflow/envs/repeatmasker.yaml` (`setup_famdb`,
`export_dfam`, `repeatmasker_chunk`, `divergence`) now put their own env
first (`ENV_PATH_GUARD` in the Snakefile) and log which `python3`,
`famdb.py` and `RepeatMasker` they resolved. Either shell state is safe for
them. `CONDA_OVERRIDE_BIN` passes through to jobs unchanged; it only
chooses which conda builds and activates the envs.

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

### Divergence landscape

`{outdir}/summary/divergence_landscape.tsv` (and its plot) is **not**
RepeatMasker's `.divsum`. `calcDivergenceFromAlign.pl` sums every
alignment in `.align`, which keeps overlapping alignments that `.out`
resolves. Across tandem arrays, adjacent monomers' alignments overlap, so
the `.divsum` can count the same bases two or more times (on the first
*M. limosa* run, the Unknown class summed to about 2× its `.out` bp).
`summarize_rm.py` instead computes a Kimura distance per alignment from its
transitions and transversions (no CpG adjustment, same as `-noCpGMod`) and
gives each genome base to the highest-scoring alignment covering it. Each
base counts once, so a class's landscape sums to about its bp in
`class_composition.tsv`; the summarize log warns past a 5% difference.
Landscape bins are 1% wide; `pct_non_n` is bp / non-N assembly length.
Simple_repeat and Low_complexity take part in the overlap resolution but are
left out of the landscape, because divergence from a consensus means nothing
for them. RepeatMasker's own `.divsum` and landscape `.html` are still
written under `{outdir}/{arm}/{species}/divergence/` for reference.

`repeatmasker.scatter_count` (default 10) and the `repeatmasker` resource
block are **per chunk now**, not per whole genome — untuned placeholders,
adjust both from real per-chunk runtimes observed in the wiring test.

## Threading summary

| step | parallelism |
|---|---|
| `BuildDatabase` | none |
| `RepeatModeler -threads` | partial; plateaus (RepeatScout/RECON serial); rounds only (no `-LTRStruct`) |
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
skipped bp, any `curated_override` in effect, and a full `config.yaml`
snapshot.
