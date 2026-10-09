# repeat_compare

RepeatModeler2 + RepeatMasker cross-sample repeat comparison, built for
*Eptatretus stoutii* (de novo) vs *Myxine limosa* (published), but
sample-agnostic and manifest-driven.

## Design

Both genomes are masked with **two arms**:

| arm | library | role |
|---|---|---|
| `shared` | non-redundant union of all samples' de novo families (+ optional Dfam export) | **primary comparison** |
| `own` | that sample's de novo families only (+ same optional Dfam export) | sanity check / concordance |

Masking every sample with only its own library biases the comparison
(each genome is best-annotated for its own families), so `shared` is the
one to trust for cross-sample numbers; `own` exists to sanity-check it via
`{outdir}/summary/arm_concordance.tsv`. (`outdir` is `results_v2` by
default: the restructured pipeline below writes to a fresh directory so
runs made with the previous `-LTRStruct` version in `results/` stay
untouched. Paths below are written as `{outdir}/...`.)

Hagfish undergo programmed germline-to-soma genome rearrangement — a
germline assembly and a somatic assembly are different genomes. This
pipeline does **not** act on the manifest's `tissue` column; it only
carries it into every summary table and warns (never fails) if samples
disagree or are `unknown`. Assembly-quality covariates (contig count,
total/N/non-N length, N50) are reported alongside every repeat result for
the same reason — a fragmented or collapsed assembly undercounts repeats.

**Naming: "sample", not "species".** Each manifest row is one assembly,
which is often one individual of a species (two meadowlark individuals,
say). The wildcard, the script arguments and the first column of every
summary table therefore say `sample` (and `sample_id`). Tables written
before this rename say `species`; to read them with newer tooling run
`sed -i '1s/^\(arm\t\)\?species/\1sample/' table.tsv`. The config keys
`summary.plot_species_order` and `library.species_prefix_sep` are now
`plot_sample_order` and `sample_prefix_sep`; the old names stop the
workflow with a message.

## Manifest

`manifest.tsv` is header-driven (manifest v2): the first non-`#` line names
the columns, in any order. Every cell needs a value with no spaces; write
`NA`, `.` or `unknown` for "not known".

| column | values | used for |
|---|---|---|
| `sample_id` | `[A-Za-z0-9]+`, unique | wildcard, library-name prefix |
| `taxon` | species name, e.g. `Sturnella_neglecta` | cross-taxon cluster counts in `discovery_summary`; satellite cross-check inference for samples without TideCluster |
| `taxid` | NCBI taxid or `NA` | provenance |
| `fasta` | path (may be gzipped) | input |
| `tissue` | `germline` / `soma` / `unknown` | carried into every summary; warning if mixed |
| `sex` | `ZW` / `ZZ` / `XX` / `XY` / `unknown` | `assembly_covariates.tsv`; warning if mixed or unknown (W/Y-linked repeats) |
| `assembly_type` | `haploid` / `primary` / `hap1` / `hap2` / `dual_hap` / `unknown` | `dual_hap` (both haplotypes in one assembly, e.g. unphased Verkko) doubles that sample's absolute copy thresholds: `family_tandem.min_copies` / `major_min_copies` / `major_min_bp`, `classify.host_max_copies` and the cross-check's eligibility bp. Percentages are unaffected |
| `accession` | token or `NA` | metadata only |

`haploid`, `primary`, `hap1` and `hap2` all mean one copy per locus and behave
the same; the distinction is provenance. A legacy 5-column manifest stops the
workflow with your rows converted to the new layout, ready to paste.

## Pipeline flow

Per sample unless noted:

```
prep_genome -> genome_fingerprint, assembly_stats
build_db -> repeatmodeler (RECON/RepeatScout rounds only, no -LTRStruct)
                     |  rounds.consensi.fa / rounds.families.stk (unclassified)
LTR side pipeline:  ltr_group_genome
    -> ltr_harvest_group + ltr_finder_group (per group) -> ltr_gather -> ltr_pipeline
merge_families (rounds + LTR, RepeatModeler's own cd-hit merge)
    -> classify_families (RepeatClassifier)
prefix_library -> cluster_library (all samples) -> shared / own libraries (+ Dfam)
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
RepeatClassifier. Every sample goes through the same path, so results are
comparable across samples. Relative to a stock `-LTRStruct` run, two things
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
- `{outdir}/{sample}/ltr/groups/manifest.tsv` shows which contigs are in
  which group;
- `{outdir}/{sample}/ltr/skipped_windows.tsv` shows which regions timed
  out;
- the per-group `*.timeouts.tsv` files list every salvaged or skipped piece.

`{outdir}/summary/ltr_skipped_composition.tsv` says what the skipped
regions hold, per sample and tool: shared-arm bp by class, by
`tandem_family` vs dispersed (any class, so an LTR-labelled satellite
counts as tandem), and the top families, each next to its genome-wide %
and the enrichment. Mostly tandem bp means the window was a satellite
array with no LTR candidates to lose; mostly dispersed LTR bp means
discovery was actually lost there.

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

The rule (`workflow/scripts/rm_extend.sh`) aims for 5 + `repeatmodeler.extra_rounds`
rounds:
- **No `RM_*` yet:** a fresh run, with `-numAddlRounds extra_rounds` when it's > 0.
- **An interrupted `RM_*`** (no `consensi.fa.classified`): it's recovered in
  place up to the target.
- **A finished `RM_*` with enough rounds:** it's reused.
- **A finished `RM_*` with fewer rounds:** it's copied to `RM_*.ext`, and the
  copy is extended with `-recoverDir`. The original is never modified, and
  rounds 1–5 of the copy must stay byte-identical to it (checked before and
  after). Raising `extra_rounds` again later extends the same `.ext`. Delete
  `RM_*.ext` to start the extension over.

**Rule of thumb:** RepeatModeler and the LTR side pipeline must both run
on the identical genome FASTA for a sample, every
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

Two tables answer "should we sample more deeply?":
- **`{outdir}/summary/discovery_round_novelty.tsv`** (the one to decide
  on). Each round's families are clustered with every earlier family
  (earlier rounds plus the LTR families) by RepeatModeler's own redundancy
  rule, cd-hit-est 80% identity over 80% of the shorter sequence. A family
  that clusters with an earlier one is a `refinement` (a sharper or longer
  consensus of something already in the library). Every other family is
  `new`. The bp each set owns comes from the own-arm `family_tandem.tsv`.
  - **Stopping rule:** another round is worth it while the last round's
    `new_dispersed_pct_non_n` is still ≥ ~0.5% of the genome.
    `refinement_pct_non_n` doesn't count: those bases were already masked
    by the earlier family. `new_tandem_pct_non_n` doesn't count either: one
    newly found satellite array says nothing about missed TE families.
- **`{outdir}/summary/discovery_round_saturation.tsv`** (kept for
  continuity): own-arm masked bp by the round that discovered each family
  (`rnd-1` … `rnd-N`, `ltr`, `other`, where `other` means Dfam and
  RepeatMasker's own simple/low-complexity calls), split into `tandem_bp`
  and `dispersed_pct_non_n`. It overstates late rounds, because a
  refinement takes over bases its earlier family already masked.
- **The gold standard, not built:** mask with and without the last round's
  families and take the difference. That's one extra RepeatMasker run per
  sample.

If sampling hasn't saturated, set `repeatmodeler.extra_rounds: 1` (the
same value for every sample) and rerun. The finished rounds are kept, and
only the new round runs (see the restart section above). Raise it one
round at a time, checking the novelty table in between. `-numAddlRounds`
in `extra_args` is refused.
- **Each extra round samples another 270 Mb** (the cap), not 3× more.
  - A ~1.2 Gb haploid bird gains about 23% of its non-N bp per round.
  - A 2.47 Gb dual-haplotype assembly (unphased Verkko) gains about 11%
    per round, and both alleles of a locus compete for the sample. Expect
    to need more rounds there before the last one stops adding families.
- **How RepeatModeler 2.0.9 resumes** (from its source; `rm_extend.sh`
  handles it):
  - `-recoverDir` needs an existing, empty `round-(h+1)/` after the last
    good round h. Without it, it prints "appears to contain a successful
    run" and stops.
  - It appends to the top-level `consensi.fa` / `families.stk`.
  - When it resumes at the 270 Mb cap, `-numAddlRounds N` runs N + 1
    rounds, so the rule passes `5 + extra_rounds − h − 1`.
  - After a resume, sampling starts from a fresh shuffle. Blocks from
    rounds 1–5 can be drawn again, but already-modelled repeats are masked
    out of each sample first.
- Each extra round costs about as much as the most expensive round, and
  later rounds mostly add low-copy `Unknown` families.
- Prefer extra rounds over a larger `-genomeSampleSizeMax`: RECON's
  all-vs-all cost grows superlinearly with sample size.

### Discovery summary

`{outdir}/summary/discovery_summary.tsv` follows the family counts from
discovery to the shared library, with one row per sample plus `ALL`.
It's built from library-stage files only, so `library_only` produces it too.

| Column | Meaning |
|---|---|
| `rounds_families`, `ltr_families` | families from RepeatModeler's rounds / the LTR side pipeline |
| `merge_removed`, `merged_families` | round families dropped as redundant with an LTR family, and what's left (`merge_families`) |
| `merged_ltr_families`, `putative_subfamilies` | LTR families in the merged set; round families tagged "putative subfamily of" |
| `classified_families`, `unknown_families` | after RepeatClassifier; those labelled `Unknown` |
| `clusters`, `sample_only_clusters`, `shared_clusters` | cross-sample cd-hit-est clusters (`cluster_library`) containing this sample's families |
| `families_in_sample_only_clusters`, `families_in_shared_clusters`, `pct_families_in_shared_clusters` | where this sample's families landed |
| `shared_clusters_label_conflict` | shared clusters whose members' `Class/Family` labels disagree (see `library_membership.tsv`) |
| `dfam_entries` | Dfam families appended to the libraries, counted by unique name (`ALL` row only) |

In the `ALL` row, family counts are summed over sample and cluster counts
are over the whole clustering. A sample's family count must equal its
members in the `.clstr`; otherwise the rule fails, since the inputs would
come from different runs.

`discovery_summary_by_class.tsv` splits families into sample-only vs
shared clusters per Class (the part of each family's own label before
`/`), for each sample and `ALL`. Sharing is usually very uneven across
classes, so read this alongside the overall percentage.

### Overlap with Dfam

RepeatClassifier compares families with the Dfam in the container, but it
only assigns a class. It doesn't record which Dfam family matched, and the
Dfam export is appended to the libraries without being compared with the
de novo families. `dfam_overlap` (run when `library.include_dfam` is set)
fills that gap. It runs `cd-hit-est-2d` of each sample's families against
the Dfam export, with `cluster_library`'s identity threshold, on both
strands:
- `{outdir}/summary/dfam_overlap.tsv`: per sample and `ALL`, and per
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

Every family in the shared-arm `.out` of each sample is checked against
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
  per family and sample, largest first: `owned_bp` (the bp it holds under
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
the shared library (`{outdir}/classify/`). These are the sample-prefixed
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
4. **host_protein**: the best host protein covers at least `min_host_cov`
   and the family has at most `host_max_copies` (50) `.out` hits in every
   sample, giving `Other/host_gene`. These are excluded from TE classes and
   counted in Other. Above that copy number the match is a TE open reading
   frame that the genome annotation called a gene ("uncharacterized LOC...",
   no TE domain, so the keyword and TEsorter filters miss it): the family
   stays Unknown with `note = host_annotated_te_orf`. On the first
   *E. stoutii* / *M. limosa* run all 65 host-protein calls had >= 100 copies
   (median ~1,400).

The family stays Unknown, with the reason in `conflict`, if:
- REXdb and GyDB disagree on order (except a REXdb Maverick/Polinton call
  against a GyDB LTR call: GyDB has no Maverick models and calls the
  Maverick integrase Gypsy-like, so REXdb's call is kept; verification
  treats it the same way), or
- a domain call coincides with a qualifying host or Rfam hit (e.g. a
  domesticated TE gene).

`diamond_best` and `rfam_best` show the best hit even below the coverage cut.
For example, a tRNA-headed SINE shows a partial tRNA hit.

`family_composition.tsv` gains `bp_from_unknown`: the bp of each
Class/Family that came from relabelled Unknown families. Report how much
moved with it, together with the `evidence` column.

**Host proteins.**
- `annotation_proteins` takes the sample's gene annotations, e.g. NCBI's
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
Families are weighted by the bp they hold in each sample (`owned_bp` from
the shared-arm `family_tandem.tsv`).

| status | meaning |
|---|---|
| `agree_superfamily` | the domain call matches the label (`LTR/ERV1` vs `LTR/ERV` counts) |
| `agree_order` | same order, different or unspecified superfamily (e.g. `LTR/Gypsy` vs Copia domains) |
| `disagree_order` | the domains point to another order, e.g. a `SINE/Alu` with Gypsy domains |
| `retroposon_vs_line` | a `Retroposon` label with LINE domains: may be an autonomous LINE, or a non-autonomous element carrying LINE fragments |
| `domain_conflict` | REXdb and GyDB disagree on order |
| `host_protein` | no domain, but the consensus matches a TE-free host protein and the family has <= `host_max_copies` hits: possible gene family or contamination |
| `host_annotated_te_orf` | no domain, matches a host protein, but > `host_max_copies` hits: a TE ORF annotated as a gene. Counts as TE support, not a disagreement |
| `rfam` / `agree_rfam` / `rfam_other_rna` | a structured-RNA hit on a non-RNA label / on a matching RNA label / on a different RNA label |
| `domain_on_rna_label` | an RNA label with TE domains |
| `no_evidence` | the screens are silent |

**`no_evidence` is not a failure.** SINEs, MITEs, solo LTRs, satellites
and decayed copies have no protein domains. Read it as "not confirmed",
never as "wrong".

Outputs:
- `{summary}/class_verification.tsv`: per sample, class and status, the
  families present, the bp they hold, and that bp as a % of the class.
- `{summary}/class_disagreements.tsv`: families with a disagreeing status
  (`disagree_order`, `retroposon_vs_line`, `domain_conflict`,
  `host_protein`, `rfam`, `rfam_other_rna`, `domain_on_rna_label`), plus
  tandem families with a TE label. Sorted by bp, with each screen's best
  hit. Review this list by hand, starting from the top.

## Library sources (`library_source.tsv`)

`{summary}/library_source.tsv` splits each sample's masked bp by where the
family holding it came from:

| source | meaning |
|---|---|
| `own_denovo` | families discovered in this sample (`<sample_id><sep>...`) |
| `denovo:<sample>` | families discovered in the other samples |
| `dfam` | the Dfam export (`library.include_dfam` / `dfam_taxon`) |
| `rm_builtin` | RepeatMasker's own simple-repeat / low-complexity screen |
| `total` | all of the above |

Columns:
- `dfam_match`: for de novo sources, `all`, then `known` (the de novo
  consensus matches a Dfam entry in `library/dfam_overlap/dfam_matches.tsv`)
  and `novel` (no match); `.` on other rows and when no matches file is used
- `n_families`: families holding any bp
- `n_hits`: `.out` hit lines (copies and fragments), the N for a methods table
- `owned_bp`
- `pct_of_masked`
- `pct_non_n`: bp over the non-N assembly length
- `bp_weighted_median_div`: how old each source's matches are

bp follow the same rule as `class_composition.tsv`: every base counts once,
for its best hit, so the `total` row matches the composition total.

The known/novel split needs `dfam_matches.tsv`. The full run (`all`) builds
it first. `report_shared_only` never builds or tracks it (so it can't pull
the discovery rules), but reads it if it is already there, so build it once
and re-run the report:

    ./runsnake 75 --configfile config.yaml --rerun-triggers mtime -- results_v3/library/dfam_overlap/dfam_matches.tsv
    ./runsnake 75 ... --forcerun library_source -- report_shared_only

`known` is a lower bound on known material: a match needs
`library.cdhit.identity` identity over `library.cdhit.coverage_short` of the de novo family (0.8 / 0.8 by default), so diverged or partial
matches count as novel, and a shared cluster counts as known only when its
representative (the name in the library) matched.

Caveats:
- Each base counts toward the family that *won* it at masking time. The Dfam
  export is appended to the library after the de novo families are
  clustered, so a Dfam entry and a near-identical de novo family both stay
  in the library and compete base by base. Where the de novo consensus fits
  better, the bp count as de novo even if Dfam has the same element.
  `dfam_overlap.tsv` lists de novo families that match known Dfam families.
- A family shared by both samples is one cd-hit cluster, kept under its
  longest member's name, so `denovo:<other samples>` includes shared
  families whose representative came from the other samples
  (`library_membership.tsv` has the clusters).
- On the first *E. stoutii* / *M. limosa* run, Dfam held under 1% of masked bp
  in both samples (old, partial matches, ~24% divergence), while about 11% of
  each genome was masked by the *other* samples' de novo families.

## Curating a family (`curation_candidates.tsv`)

`{summary}/curation_candidates.tsv` ranks families worth a manual recheck:
those holding >= `summary.curation_min_pct_masked` (0.1%) of a sample's
masked bp and flagged as
- `unknown` / `unknown_tandem`: still Unknown after the classify screens
- `long_for_class(len>max)`: consensus longer than plausible for its label
  (e.g. > 1.5 kb for a SINE, > 5 kb for Tc1/hAT/PiggyBac, > 9 kb for a LINE),
  usually a chimeric consensus or a wrong label
- `verify:<status>`: a class_disagreements.tsv status (Maverick
  `domain_conflict` is expected and not flagged)
- `young(div%)`: median `.out` divergence < `summary.curation_young_div`
  (3%), a recent burst whose near-identical copies are easy to curate
- `ltr_pipeline_dna_label`: found by the LTR pipeline (`<sp>_ltr-N_family-M`,
  a structural LTR candidate) but labelled `DNA/...`; the label is probably
  wrong

Families already in `summary.element_groups` or `classify.curated_families`
are listed last with their element. The recheck that resolved `Esto_rnd-1_family-332` (a non-autonomous
LTR element split across 26 families):
1. `workflow/scripts/family_profile.py` (consensus structure, copy
   coverage, full-length copies, TSDs), then again with `--flank 5000`
   to see whether copies continue past the consensus ends
2. align the extended copies and build a consensus (e.g. abPOA, or MAFFT
   plus a gap-aware majority consensus)
3. `family_groups.py members` to list the library families that are pieces
   of it; `ltr_ages.py` for LTR elements
4. name it (below) and add its rows to `classify.curated_families`

## Naming curated elements

Curated elements follow the Dfam / RepBase convention:

    <Superfamily>-<n>[N]_<sample>      e.g. Gypsy-1_Esto, Gypsy-N1_Esto, hAT-N2_Esto

- `<Superfamily>`: the RepeatMasker superfamily (Gypsy, Copia, ERV1, CR1,
  RTE, hAT, TcMar, PiggyBac, Maverick, Helitron, ...), or the class when
  the superfamily is not established (`LTR-N1_Esto`)
- `<n>`: running number per superfamily and sample, in the order elements
  are curated; never reused
- `N`: non-autonomous (no coding capacity of its own; mobilised in trans),
  as in RepBase's `hAT-N1_DR`
- `_<sample>`: the pipeline's sample ID (`_Esto`, `_Mlim`), matching
  Dfam's `Naiad_Ebur` style
- LTR elements: the curated library entries are `<name>-LTR` and `<name>-I`
  (internal region), as RepeatMasker libraries store them; the table's
  `part` column records which library families are which
- satellites: `SAT-<n>_<sample>` (Dfam style), with the monomer length in
  `note`

The curated-families table maps the original library families (kept as the
permanent IDs for traceability) to the element and its RepeatMasker
`class_family` (`LTR/Gypsy`, never the element name): the format of
`family_groups.py members` output, with `class_family` filled
(`members --class-family LTR/Gypsy --ltr-len 1695`). Rows from several
elements go in one file; when elements are related (a non-autonomous
element and its partner share LTRs and other pieces), run `members` per
element and `family_groups.py resolve` to keep each library family under
the element it matches best. Curated consensi live in `curation/elements/`.
Set `classify.curated_families` to it:
- `summarize_rm.py` labels every member family's hits with `class_family`,
  over RepeatClassifier, the tandem carve-out and the classify screens
  (`family_composition.tsv` column `bp_curated` shows how much)
- `element_groups.tsv` sums each element's members
- `curation_candidates.tsv` lists them as curated
Nothing is remasked; rerun the report with `--forcerun summarize` after
editing the table.

## Element groups (`element_groups.tsv`, optional)

RepeatModeler often splits one element across several library families,
e.g. an LTR retrotransposon whose internal region is one `rnd` family and
whose LTRs are several `ltr` families. To report the element's whole
footprint:

1. Build a curated full-length consensus (e.g. extend the family with
   `workflow/scripts/family_profile.py --flank 5000`, align the copies,
   take a majority consensus).
2. `blastn -query element.fa -db <library> -outfmt "6 qseqid sseqid pident length qstart qend sstart send qlen slen evalue bitscore" > element_vs_lib.tsv`
3. `python3 workflow/scripts/family_groups.py members --blast element_vs_lib.tsv --group <name> --out groups.tsv`
   lists library families with >= 80% identity over >= 50% of their own
   length. Review the list; concatenate several elements' outputs into one
   file.
4. Set `summary.element_groups: groups.tsv` and rerun the report. (For a
   named element with a curated class, use `classify.curated_families`
   instead: it also relabels, and its elements are grouped the same way.)

`element_groups.tsv` has one `group` row per sample and element (members
present, `.out` hits, bp, % of masked, % non-N) and one `member` row per
family. Each base still counts for the family that won it, so a group is the
sum of its members' `owned_bp`, and nothing is remasked.

Insertion ages of an LTR element come from the LTR-pipeline candidates:
`workflow/scripts/ltr_ages.py` takes the pipeline's `{sample}/ltr/rawLTR.scn`,
keeps candidates whose internal region is covered by the element's
families, and reports the 5'/3' LTR similarity of each intact copy (with
`--rate`, ages via the Jukes-Cantor distance, T = d / 2r).

## Placing Unknown families by their neighbours (`family_neighbors`)

Many Unknown families are pieces of elements whose other pieces are
classified: the LTRs of a Gypsy internal region, the 5' end of a LINE, a
split consensus. Their copies sit next to those pieces again and again, on a
fixed side and in the same orientation. `family_neighbors.py` counts each
family's neighbours on its 5' and 3' sides (in the family's own orientation)
in the shared-arm `.out` files and calls the patterns below. Neighbours are
looked for in widening windows (`summary.neighbors.max_gap`, default 100,
500 and 2000 bp between the copy and its neighbour). Each family takes its
call from the tightest window that gives one; pieces of young elements
abut, while old copies are split by unmasked, diverged stretches. If a
wider window finds the same pattern and partner with more support, that
support is reported, with its window in `max_gap`. Because a 2 kb window in
a repeat-dense genome nearly always holds some neighbour, every pattern
also needs the partner on the same strand ≥ `strand_ratio` times as often
as on the opposite one (random neighbours are 50/50), and calls found only
in windows wider than `high_max_gap` are capped at medium.

| pattern | copies look like | proposal |
|---|---|---|
| `LTR_of` P | P (same strand) on the family's 3' side in some copies, 5' side in others; the family also sits at both ends of P's copies | P's class, `part = LTR` |
| `internal_of` P | P on both sides, same strand | P's class, `part = I` |
| `5prime_of` / `3prime_of` P | P continues the family on one side, same strand | P's class, `part = 5prime` / `3prime` |
| partner-anchored | as above, but seen from P: ≥ `anchor_frac` of P's copies carry the family on the facing side, though the (much larger) family is mostly elsewhere. P must be one of this library's own families (sample-prefixed) with ≥ `anchor_min_copies` copies: a foreign Dfam consensus with a few dozen hits is usually a partial match to the family itself | P's class, capped at medium; used only when the family's own copies show no pattern |
| `internal_of_LTRs` / `next_to_LTRs` | LTR variants pooled: LTR-pipeline (or LTR-labelled) families, summed, sit on both sides (or one side) of the family, same strand. Catches internal regions whose LTRs are split across many variant families, none frequent enough alone | `LTR/Unknown` (or the LTR class carrying most of the pooled support), `part = I`; one side only: medium |

- Targets are families labelled `Unknown` or `<Order>/Unknown` after the
  curated, cross-check and reclassification labels, with ≥ `min_copies`
  copies. A call needs ≥ `min_pairs` copies and ≥ `min_frac` of them;
  `high` needs ≥ `high_frac` and twice `min_pairs`.
- An `<Order>/Unknown` partner (e.g. LTR/Unknown LTRs) lifts a plain Unknown
  to that order. A target linked only to another Unknown gets no class but a
  note; if that partner is called in the first pass, the target inherits
  the call (medium).
- Orientation: RepeatModeler and LTR-pipeline consensi are oriented
  arbitrarily, so one piece of an element can be stored reverse-complemented
  relative to another. A partner consistently on the opposite strand counts
  like one on the same strand (the note says "reverse strand"); what is
  required is that one orientation dominates (`strand_ratio`). Pooled LTR
  partners each count in their own dominant orientation, and only if that
  orientation is significant (binomial p ≤ `pool_max_p` against 50/50) and
  the partner is substantial (≥ `pool_min_copies` copies and ≥
  `pool_min_share` of the family's copies): in LTR-rich genomes many small
  partners passing a ratio test by chance would otherwise add up to a call.
- Only TE orders (DNA, LINE, SINE, LTR, RC, Retroposon, PLE) and rRNA
  (rDNA spacer pieces) pass their label on; `RNA`, `tRNA`, `Other` ... and
  the labels in `no_transfer` (default `SINE/Alu`, implausible in hagfish)
  give a note instead.
- Unknown families that keep joining each other (neither side has a usable
  label) are grouped into chains in `family_neighbors_chains.tsv`: likely
  one element split across several consensi, to profile or extend as one.
- Notes only: `3' poly(A)` (an A-rich simple repeat right after the 3' end:
  non-LTR retrotransposon or SINE) and `tandem` (copies next to copies of
  the same family).
- Random insertions near other repeats give mixed partners and 50/50
  strands, so they rarely reach `min_frac`.

Outputs in `summary*/`: `family_neighbors.tsv` (top partners per side),
`family_neighbors_calls.tsv` (one row per tested family) and
`family_neighbors_chains.tsv`,
`family_neighbors_proposals.tsv` (`classify.curated_families` format;
`group` is the partner's curated element name when it has one, else
`with-<partner>`). Report-only: copy accepted rows into your curated table.
Built by `all` and `report_shared_only` (`summary.family_neighbors: false`
turns it off).

## TEtrimmer (`tetrimmer_only`)

[TEtrimmer](https://github.com/qjiangzhao/TEtrimmer) automates the manual
curation steps (copy extraction, extension, MSA cleaning, boundary,
LTR/TIR/poly(A)/TSD and Pfam checks, RepeatClassifier) and writes a PDF
report per consensus. Here it runs only on the largest families still
unresolved, to classify them, not to replace the library:

1. `tetrimmer_select`, per sample: families still `Unknown` /
   `<Order>/Unknown` after the curated, cross-check and reclassification
   labels and the `family_neighbors` calls at `tetrimmer.neighbor_confidence`,
   not tandem in any sample, whose largest footprint is in this sample: the
   top `tetrimmer.top_n` by owned bp with ≥ `tetrimmer.min_bp`. Each family
   runs once, on the genome where it is largest.
2. `tetrimmer`, per sample, conda env `workflow/envs/tetrimmer.yaml`
   (`--classify_unknown` by default; `tetrimmer.preset`, `classify`,
   `pfam_dir`, `extra_args`). Resources: `resources.tetrimmer` (32 threads,
   120 h by default; TEtrimmer's own benchmark was ~2-5 h for 3,500-8,600
   families on a 1.7 Gb genome with 48 cores, and hagfish families have far
   more copies). The run directory `tetrimmer/{sample}/run` is not a
   declared output: a killed job resumes with `--continue_analysis`; delete
   it to start over.
3. `tetrimmer_report`: `tetrimmer/tetrimmer_summary.tsv` (every output
   consensus, mapped back to library family IDs) and
   `tetrimmer/tetrimmer_proposals.tsv` (curated-table format, one row per
   family whose label TEtrimmer changed, from its best output: Perfect →
   high, Good → medium, else low; split families noted).

Review before accepting: the PDFs in
`tetrimmer/{sample}/run/TEtrimmer_for_proof_curation/` (or TEtrimmerGUI, which
needs only Python), then copy rows into your curated table.

    ./runsnake 4 --configfile config.yaml --config mask_shared_library=results_v3/library/shared_library.fa \
        --rerun-triggers mtime -- tetrimmer_only

Notes:
- Pfam: TEtrimmer downloads it on first use. Compute nodes without
  internet need `tetrimmer.pfam_dir` pointing at a local copy.
- RepeatClassifier inside the TEtrimmer env uses that env's RepeatMasker
  libraries (the Dfam root partition bioconda ships), not the pipeline's
  FamDB setup.
- The env is unpinned (`tetrimmer` from bioconda); after the first build,
  pin the resolved version in `workflow/envs/tetrimmer.yaml`. Upstream notes
  the bioconda package can lag GitHub.

## Satellite, rDNA and mito cross-check (`satellite_crosscheck`)

`family_tandem` infers arrays from RepeatMasker hits alone. This
cross-check tests library families against three independent tools:
- TideCluster (tandem-repeat arrays clustered into TRCs);
- ribotin (rDNA unit models);
- MitoHiFi (mitogenomes);
- optionally, the satellite pipeline's harmonized motif library
  (`harmonized_library`, e.g. `harmonized_repeatmasker_lib.fasta`).

The two satellite sources answer different questions, and either works
alone. TideCluster says whether a family *sits in arrays in this assembly*;
the harmonized library says *which known satellite unit* a consensus is
made of, and names it. Without the harmonized library (before the full
satellite pipeline has run), calls rest on TideCluster alone, as before.
The library also holds TE-derived and dispersed "satellite" units (the
reason the old satellite arm was removed, below), so a motif match alone
never makes a satellite: see the motif rule under **Calls**.

It is **report-only**: it writes evidence and proposals, and nothing is
relabelled or merged.

**Inputs.** Point `satellite_crosscheck.external_annotations` at a TSV with
the header `sample_id tidecluster_dir tidecluster_prefix ribotin_fa
mitohifi_fa` (see `external_annotations.example.tsv`).
- `sample_id` is the manifest's `sample_id`.
- Every cell needs a value. `.`, `NA`, `na`, `no`, `false` and `none` mean
  "not provided"; blanks and spaces are refused.
- From `{tidecluster_dir}/{tidecluster_prefix}_*` (TideCluster 1.21.3
  `run_all`), the pipeline reads:
  - required: `clustering.gff3`, `seqid_lengths.tsv`,
    `consensus_dimer_library.fasta`;
  - if present: `tidehunter.gff3` / `tidehunter_short.gff3` (per-array
    monomers), `tarean_report.tsv` (TRC monomer length),
    `trc_superfamilies.csv` and `rdna.tsv`.
- **Same assembly required.** TideCluster must have run on the same
  assembly as the manifest FASTA. `trc_regions` compares
  `seqid_lengths.tsv` with the pipeline's `fingerprint.tsv` (through
  `name_map.tsv`) and stops on any mismatch. To skip the check, set
  `allow_seqid_mismatch: true`; the mismatched contigs are then dropped.

**Run it on finished results** without redoing anything upstream:

    snakemake satellite_crosscheck_only --rerun-triggers mtime -n   # should list only the new jobs
    snakemake satellite_crosscheck_only --rerun-triggers mtime ...  # usual runsnake arguments

With `external_annotations` set, `all` builds it too.

**Outputs:**

| file | content |
|---|---|
| `shared/{sample}/satellite/trc_crosscheck.trc.tsv` | per TRC: arrays, bp, bp RepeatMasker masks / misses, the families holding it. Row `ALL`: total array bp the library misses, and TideHunter arrays outside any TRC |
| `shared/{sample}/satellite/trc_crosscheck.family.tsv` | per family: share of its bp inside TRC arrays, which TRCs, median TideHunter monomer length of the arrays it sits in |
| `summary/satellite_family_evidence.tsv` | family × sample. Samples without TideCluster listed in `infer_from` get the donor's TRC evidence, marked `inferred:<donor>` |
| `summary/satellite_family_calls.tsv` | one row per family: proposed class, confidence, reason |
| `summary/satellite_family_pairs.tsv` | pairs of tandem families: consensus similarity, shared-TRC fraction, monomers, verdict |
| `summary/satellite_proposals.tsv` | `classify.curated_families` format; high and medium rows only |
| `summary/satellite_trc_pairs.tsv` | TideCluster consensi matching across samples: coverage both ways, identity, KITE founder of each. These are the satellites shared between samples |

**Calls**, in priority order:
1. `rRNA` (high): ≥ `rdna_min_cov` of the consensus matches a ribotin model
   at ≥ `rdna_min_id` %.
2. `Other/NUMT` (high): the same test against the mitogenome, with the
   `mito_*` thresholds.
3. `Satellite`:
   - **high:** ≥ `min_trc_cov` of its bp in TRC arrays (that alone is
     array evidence; the `tandem_family` flag is not needed, since long
     units with few copies per array can miss its thresholds), in every
     *eligible* TideCluster sample. A sample is eligible when the family
     holds ≥ `family_tandem.major_min_bp` there AND is tandem there or
     holds ≥ `eligible_min_share` of its bp in its top sample, so a few
     scattered copies in the other species don't veto the call;
   - **medium:** only some of those samples pass;
   - **medium (partial):** tandem, with between `partial_trc_cov` and
     `min_trc_cov` of its bp in TRCs. Typically these are arrays only partly
     clustered by TideCluster.
   - **low:** TRC support only in samples below that size.
   - **Two different period limits in TideCluster:**
     - TideHunter's `-P` limits which arrays are *found*: 3000 by default
       (TideCluster `-T`), 25000 with `--long`.
     - KITE (`{prefix}_kite/monomer_size_top3_estimats.csv`) then re-measures
       each found array's own *founder period*, up to 10 kb (25 kb when
       extended). That is the x-axis of TideCluster's cluster-overview plot,
       and why founders > 3 kb appear. Such arrays were found through a
       sub-period ≤ `-P`, e.g. a HOR's basic monomer.
     - A unit with no internal period ≤ `-P` is never found.
     - KITE founders are read into `trc_info.tsv`
       (`kite_founder_median`, weighted by array length) and per family
       (`kite_founder_median`). The pair `monomer_match` uses them first,
       then the TideHunter median, then the consensus self-period. A
       "multiple" must be within `monomer_tol` of the shorter unit.
   - **TideHunter's detection limit:** TideHunter only looks for repeat units up to
     `-P` bp, which is 3000 by default (TideCluster `-T`) or 25000 with
     `--long`. The limit is read from `cmd_args.json` into
     `shared/{sample}/satellite/tidecluster_params.tsv`. A tandem family
     whose unit (`monomer_period`, else the consensus length) is longer than
     that cannot be in a TRC, so in that sample it neither passes nor fails;
     it is listed under `beyond_tidehunter`. If its `tandem_frac` is
     ≥ `min_tandem_frac_long`, it gets `Satellite` on the RepeatMasker arrays alone:
     medium for `Unknown` families, low for TE-labelled ones, since those may
     be tandem segmental copies of the TE. The unit is the KITE founder
     when there is one. Rerunning TideCluster with `--long` turns
     these into a real TRC test.
   - A TRC that TideCluster flags as rDNA turns a Satellite call into
     `rRNA` (medium).
4. **Harmonized motifs** (with `harmonized_library`). Each motif is tiled
   (rotation-safe) and the library consensi are searched against it; the
   best motif per family is reported (`best_motif`, `motif_cov`,
   `motif_id`, `motif_tier`):
   - **high:** ≥ `motif_high_cov` of the consensus at ≥ `motif_high_id` %;
   - **medium:** ≥ `motif_min_cov` at ≥ `motif_min_id` %;
   - **partial:** ≥ `motif_partial_cov` at ≥ `motif_high_id` %: a satellite
     segment plus other sequence (chimeric consensus?). Adds `low` at most,
     and `part = partial` in the proposals.

   A high/medium motif tier can raise a Satellite call (or an empty one) to
   that tier, but only with array evidence in some sample (`tandem_family`,
   `tandem_frac` ≥ `motif_array_tandem_frac`, or `frac_in_trc` ≥
   `partial_trc_cov`). Without it the call stays `low` ("TE-derived
   motif?"); TE-labelled families with arrays stop at `medium` ("TE-derived
   satellite?"). Families RepeatMasker already labels `Satellite` skip the
   array test: the motif only names them. Proposals take the motif ID as
   the `group` name (a partial match too, with `part = partial`), so
   curated satellites carry the satellite pipeline's names (README "Naming
   curated elements"). The best motif per family maximizes coverage ×
   (identity − `motif_min_id` + 1): a near-identical long unit over part of
   the consensus beats a short motif loosely tiled over more of it.

A family currently under a TE label keeps that label in its `reason` / `note`
(e.g. `LTR/ERVK-derived`). Inferred samples can confirm a call (`ind8:confirms`)
but never make one.

**Pair verdicts.** Pairs are tested from three sources: all pairs in
`focus_families`, families sharing a TRC, and families whose consensi hit
each other. The verdicts:
- `same_satellite` needs all three of these:
  - consensus identity ≥ `min_pair_id` over ≥ `min_pair_cov`;
  - TRC sharing ≥ `min_shared_trc_frac`. The overlap coefficient is
    Σ min(a_t, b_t) / min(Σa, Σb) over TRCs t;
    only samples where both families hold ≥ `min_shared_trc_bp` in TRCs
    count. Otherwise a sprinkle of hits inside another satellite's arrays
    scores 1.00; such samples show as `low_bp(a/b)`;
  - monomers within `monomer_tol`, or an integer multiple (HOR-like).
- `same_satellite_diverged`: the same arrays and monomer, but the consensi
  only reach `related_min_id` (not `min_pair_id`). These are likely
  variants or subfamilies of one satellite. They are reported, not grouped.
- `co_located_related`: shared arrays and related consensi, but different
  monomers, so possibly a composite array.
- `same_arrays_same_monomer`: shared arrays and the same founder, but
  consensi that blastn can't align. Probably one satellite with very
  diverged consensi; check it by hand. These are never grouped.
- `co_located_distinct`: shared arrays, unrelated consensi.
- `related_not_co_located`: related consensi (≥ `related_min_id`) in
  different arrays.
- `related_no_shared_sample`: related consensi, but no TideCluster sample
  holds both.
- `independent`: consensi overlap < `independent_max_cov`, and TRC sharing
  ≤ `max_independent_shared_frac`.
- `undetermined`: anything else.

Pairs are only formed among medium/high calls and `focus_families`.

How monomers and consensi are compared:
- **Monomers:** a pair's monomers are the two KITE founders in the sample
  where the families share the most arrays (`monomer_a/b` names it). Only
  when no sample has both does each family fall back to its own best
  estimate.
- **`top_trc_match`:** lists the families' top TRCs that are the same TRC
  (`=same`) or whose consensi match across samples (`~`). Such a
  cross-sample match with agreeing founders upgrades
  `same_arrays_same_monomer` to `same_satellite_diverged`.
- **Rotation:** consensus comparisons are rotation-agnostic.
  - TideCluster consensi are dimers, and most RepeatModeler satellite
    consensi are multimers, so every rotation lies in one piece inside them.
  - blastn searches both strands, and coverage merges all HSPs.
  - A library consensus shorter than 1.5× its unit is searched as a dimer,
    with coordinates folded back to its own length.
- **Simple repeats:** RepeatMasker's built-in `Simple_repeat` /
  `Low_complexity` entries are never called.

**Applying the labels (`satellite_crosscheck.apply: true`).**
`satellite_apply` keeps the proposals at `apply_confidence` (default
`[high]`) and drops any family your `classify.curated_families` already
lists, giving `summary/satellite_applied.tsv`. `summarize` then reads both
tables, yours first, and the first table that lists a family wins.
`family_composition.tsv` shows the effect: `bp_curated` holds what your
table relabelled, `bp_crosscheck` what the applied proposals relabelled.
Nothing is remasked. Grouping (shared `group` names) never changes labels.
With `apply: false` the proposals stay report-only, and you curate by hand.

A pair listed in `expected_independent` is never grouped. If the test calls
it `same_satellite`, it is flagged `CONTRADICTS expected_independent` for a
manual look. Only `same_satellite` pairs share a `group` name in the
proposals file. Merging is never done automatically.

## Genome-wide NUMT calls (`numt_only`)

The cross-check above labels a library *family* `Other/NUMT`. Most NUMTs are
single-copy, though, and never become RepeatModeler families, so they are
not reported there. This step calls NUMTs per locus. It is **report-only**:
nothing is relabelled or remasked. A later base-level reconciliation will
give NUMT bp priority over TE labels.

**Inputs.** The `mitohifi_fa` column of `satellite_crosscheck.external_annotations`
(MitoHiFi `final_mitogenome.fasta`). Samples without one are skipped, with a
note at startup and a `no_mitogenome` row in the summary. The step also needs
the prepped genome and its `fingerprint.tsv` / `assembly_stats.tsv`, and the
shared-arm `.out`.

**Run it on finished results:**

    snakemake numt_only --rerun-triggers mtime -n   # should list only numt_* jobs
    snakemake numt_only --rerun-triggers mtime ...  # usual runsnake arguments

`all` builds it too when any sample has a `mitohifi_fa`.

**Method** (`numt_blast`, then `numt_calls`, per sample; `workflow/scripts/numt_calls.py`):
1. The mitogenome is doubled (seq+seq), so an alignment across the circular
   origin is one HSP. Hit coordinates are folded back modulo the mitogenome
   length; `mito_wraps = yes` marks a call that crosses the origin.
2. `blastn -task dc-megablast`, mitogenome as query against the nuclear
   genome, with `-evalue numt.evalue` and `-dust no`. `-max_target_seqs` is
   set above the assembly's contig count (BLAST's default of 500 would
   silently drop contigs), and there is no `-max_hsps`, so every NUMT on a
   chromosome is kept. blastn rather than minimap2: old NUMTs are short and
   70–85 % identical, below what minimap2's seeds reliably find, and
   BLAST/LAST is the usual choice in NUMT studies.
3. **Mitochondrial contigs** are excluded. A contig is mitochondrial when hits
   at ≥ `mito_contig_id` % cover ≥ `mito_contig_cov` of it. These contigs are
   listed in the log and `mito_contigs.tsv` and get no calls.
4. **Doubled-query duplicates.** The same genome locus is reported once per
   query copy, whole or as pieces at the copy edges. Taking hits by
   decreasing bitscore, a hit whose genome interval overlaps an already kept
   hit (same contig and strand) by ≥ 50 % of its own length is dropped.
5. Two levels of call:
   - **hit-level** (`numt_hits.bed`): one per remaining HSP. This is the
     primary, unambiguous count.
   - **compound** (`numts.bed`): neighbouring hits on the same contig and
     strand, joined when the *mitogenome* coordinates continue (the next hit
     starts within `mito_merge_gap` bp of where the previous one ended, as a
     gap or an overlap, across the origin too) and the *genome* gap is
     ≤ `merge_gap`. That is one NUMT split by an indel or a later insertion,
     not two insertions. Single hits are compound calls with `n_hits = 1`.
     `max_gap` is the largest genome gap. `gap_fill` is what fills it in the
     shared-arm `.out`, as `class:bp` with each base given to its best hit and
     the labels `class_composition.tsv` uses; `unmasked` means no hit.
6. Both levels keep calls with aligned bp ≥ `min_len`. The filter runs after
   merging, so a short piece can still join a compound call.

**Choosing `merge_gap`.** The default (500 bp) only bridges indels and short
insertions. `numt_gap_hist.tsv` bins the genome gap of *every* mito-colinear
neighbour pair, whatever its size, with how many gaps one class fills
(≥ 80 %) and the fill's class bp. If the gaps cluster at TE lengths (e.g.
1–6 kb, filled by one TE), raise `merge_gap` to cover that cluster and say so
in the methods. Without that evidence, don't start from 2 kb. Only
`numt_calls` reruns when `merge_gap` changes; the blastn is kept.

**Outputs.** Coordinates are on the prepped genome, with the names used in the
`.out` and every other output. BED files are 0-based half-open with a `#`
header.

| file | content |
|---|---|
| `numt/{sample}/numt_hits.bed` | hit-level calls: contig, start, end, id, aligned_bp, strand, identity, mito_start, mito_end, mito_wraps |
| `numt/{sample}/numts.bed` | compound calls: the same columns (identity bp-weighted; mito_start/end along the NUMT) + n_hits, max_gap, gap_fill |
| `numt/{sample}/numt_gap_hist.tsv` | genome gaps between mito-colinear neighbouring hits, binned (for choosing `merge_gap`) |
| `numt/{sample}/mito_contigs.tsv` | contigs excluded as mitochondrial, with cov and identity |
| `numt/{sample}/mito_vs_genome.blastn.tsv` | the raw blastn table (query coordinates on the doubled mitogenome) |
| `summary/numt_summary.tsv` | per sample: mito contigs, n hit-level / compound NUMTs, `numt_bp` (union of hit-level calls) and its % of non-N bp, `numt_span_bp` (union of compound spans), median identity and aligned-length quantiles of compound calls, `bp_overlap_<class>` (shared-arm `.out` bp inside hit-level calls, by class) |

**Sanity check.** When the satellite cross-check is configured, the families
it calls `Other/NUMT` should sit inside called NUMTs. For each one, the log
lists its `.out` bp in the sample (off mitochondrial contigs) and how much of
that falls inside hit-level calls. A family mostly outside calls is flagged.
The totals are the summary's `crosscheck_numt_*` columns (`NA` without the
cross-check).

Unit tests for the origin folding and the merge logic:
`python3 -m unittest discover tests`.

## Structural TIR / Helitron pilot (`dna_te_pilot_only`)

RepeatModeler's RECON/RepeatScout rounds model TIR elements poorly, and
Helitrons hardly at all. The LTR side pipeline already adds structure-based
LTR families; this is the DNA-TE analogue, starting as a **pilot**. The pilot
is report-only and changes nothing in the library. Adding families would
change the library, so both masking arms and everything downstream would
rerun. That integration (Stage B2) is only worth building if the pilot shows
real new bp.

    snakemake dna_te_pilot_only --rerun-triggers mtime -n   # only dna_te_* jobs on finished results
    snakemake dna_te_pilot_only --rerun-triggers mtime ...  # usual runsnake arguments, + --use-singularity

It is not part of `all`.

**Callers.** EDTA's structural modules run standalone: `EDTA_raw.pl --type tir`
(TIR-Learner) and `--type helitron` (HelitronScanner). Each module's own
filters (flank-repeat test, tandem and TEsorter cleanup) are applied, so the
output is EDTA's `*.intact.raw.fa`. They run in the EDTA image `dna_te.container`
(biocontainers `edta:2.3.0--hdfd78af_0`, pinned by digest). TIR-Learner's
model is `dna_te.tir_species` (`others`; EDTA only has rice and maize models
besides it).

**Parallelism and stalls**, mirroring the LTR tools:
- `dna_te_group_genome` writes bp-balanced whole-scaffold groups like
  `ltr_group_genome`'s, `dna_te.n_groups` of them (one job per group and type).
  They are written separately: the LTR groups are `temp()`, and re-making them
  would make every finished `ltr_*` job look out of date. Changing
  `n_groups` reruns every DNA-TE group but no LTR job.
- Each window's start and end are logged with a timestamp, its runtime and
  its candidate count. On *E. stoutii*, TIR-Learner took ~7 min per 5 Mb
  window (4 at a time, 16 threads; 2 h for a ~355 Mb group, no timeouts).
- `dna_te_candidates_group` (one job per group and type) cuts the group into
  windows of `window_size` bp. Long scaffolds overlap by `overlap`; small
  scaffolds are packed together. It runs `threads / threads_per_window`
  windows at once, each under `window_timeout_s`.
- A window that times out is killed and logged in
  `{sample}/dna_te/skipped_windows.tsv` (tool, contig, start, end, bp,
  reason), like `ltr/skipped_windows.tsv`.
- A window that fails (EDTA exits without its result) is rerun one scaffold
  piece at a time, and only the pieces that fail again are logged, with
  reason `failed exit N (tandem array | ordinary sequence, zlib R; log PATH)`.
  The class comes from the piece's zlib compression ratio: below
  `dna_te.tandem_zlib` (0.1) it is a tandem array (ordinary DNA compresses to
  ~0.25-0.30). Example: a 2.08 Mb meadowlark scaffold that is one GGAA-rich
  satellite array (zlib ~0.03) breaks TIR-Learner's TIRvish parsing; it holds
  no TIRs to find.
- Tandem-array failures are expected (bird assemblies carry large satellite
  scaffolds) and never stop the job. Failures on ordinary sequence point to a
  tool problem: by default (`dna_te.max_failed_frac: 0`) any of them fails the
  job, so it isn't hidden behind skip rows. The one seen so far (swifter's
  process scheduler inside TIR-Learner's worker pool, below) struck only the
  windows with the most TIRvish hits -- skipping those would bias the pilot.
- TIR-Learner 3 in the EDTA 2.3.0 image needs two shims, loaded into EDTA's
  Python processes only (`workflow/scripts/compat/edta_python/sitecustomize.py`):
  pandas 3 dropped the positional `row[0]` fallback `check_TIR_TSD.py` relies
  on (`KeyError: 0`), and swifter 1.4.0 parallelises a large apply with
  dask's process scheduler inside a `multiprocessing.Pool` worker
  ("daemonic processes are not allowed to have children"); it is switched to
  threads. The job log's first line says the shims are on.
  `skipped_bp` in the pilot table counts timeouts and failures per sample and
  type; `skipped_windows.tsv` says which was which.
- Whether TIR-Learner or HelitronScanner stall on the hagfish assemblies is
  still to be seen on the first run. Check the skipped windows, and tune the
  timeout from the per-window times in the log.

**Pilot** (`dna_te_pilot`, per sample and type):
1. Candidates (window-overlap duplicates dropped, `candidates.tsv` with
   scaffold coordinates) are clustered with cd-hit-est (`library.cdhit`
   settings, both strands).
2. Each cluster representative is matched with cd-hit-est-2d against the
   sample's own `-families.fa` and against the shared library.
3. `pilot_clusters.tsv` lists every cluster with the family it matches and
   that family's Class/Family (the library header label, i.e.
   RepeatClassifier's).
4. **Masked-bp proxy:** the unmatched representatives are blastn'd
   (dc-megablast) against the genome. `unmatched_cover_bp` is the union of
   their HSPs; `unmatched_new_bp` is the part the shared-arm `.out` does not
   mask yet. This is a lower bound: blastn misses old copies that
   RepeatMasker -s would find. It was chosen over a RepeatMasker run for cost.

`summary/dna_te_pilot.tsv` has one row per sample × type: candidates,
clusters, matched (own / shared / either), unmatched, the bp above, the
Class/Family counts of the matched clusters, `absent_superfamilies`, and a
`decision`. The decision rule is on the table's first line:
- integrate if the unmatched clusters add ≥ `dna_te.min_new_pct_non_n` %
  (0.5, the same scale as the round-novelty stopping rule) of non-N bp, or
- if they carry a superfamily with no family in the shared library (EDTA's
  DTA/DTC/DTH/DTM/DTT/Helitron, mapped to hAT/CMC/PIF-Harbinger/MULE/TcMar/Helitron).

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
RepeatMasker with an already-built shared library on every sample, with no
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
Output goes to `<outdir>/shared/<sample>/repeatmasker/`. Leave
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
add `--forcerun summarize`. That reruns the per-sample summaries and
everything downstream, but not the masking.

The target refuses to run without `mask_shared_library`. Without it the
shared arm would mask with `<outdir>/library/shared_library.fa`, which
schedules the whole discovery chain and re-runs RepeatMasker.

It schedules `assembly_stats`, `divergence` and `summarize` per sample,
then `combine_summaries_shared_only` and `plot_shared_only`:
- `<outdir>/summary_shared_only/`: `class_composition.tsv`,
  `family_composition.tsv`, `divergence_landscape.tsv`,
  `assembly_covariates.tsv` (same formats as `summary/`);
- `<outdir>/plots_shared_only/`: `class_composition_shared.png`,
  `divergence_landscape.png`;
- `<outdir>/shared/<sample>/divergence/<sample>.landscape.html`:
  RepeatMasker's own landscape page.

These go to separate directories so a later full run's `summary/` and
`plots/` are never mixed with them. Not produced: `arm_concordance` (needs
the own arm), `discovery_round_saturation` / `discovery_round_novelty` (own arm +
RepeatModeler rounds), `ltr_discovery`, `discovery_summary`, `dfam_overlap` and
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
`{outdir}/{sample}/repeatmodeler/{sample}.repeatmodeler_provenance.txt`
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

Each sample's genome is masked via a 3-rule scatter/gather
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

`gather_repeatmasker` merges three file types per (arm, sample):
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
written under `{outdir}/{arm}/{sample}/divergence/` for reference.

`repeatmasker.scatter_count` (default 10) and the `repeatmasker` resource
block are **per chunk now**, not per whole genome — untuned placeholders,
adjust both from real per-chunk runtimes observed in the wiring test.

## Threading summary

| step | parallelism |
|---|---|
| `BuildDatabase` | none |
| `RepeatModeler -threads` | partial; plateaus (RepeatScout/RECON serial); rounds only (no `-LTRStruct`) |
| `ltr_harvest_group` / `ltr_finder_group` | yes: `n_groups` SGE jobs per sample × `threads` windows each (threads must be ≥ 2) |
| `ltr_pipeline` (LTR_retriever, MAFFT, NINJA) | yes, `-threads`; once per sample (needs the whole genome's candidates) |
| `merge_families` (cd-hit-est), `classify_families` (RepeatClassifier) | yes, `-T` / `-threads` |
| `cd-hit-est -T` | yes |
| `RepeatMasker -pa` | yes; each slot ~4 cores under RMBlast (`repeatmasker.cores_per_pa`) |
| `calcDivergenceFromAlign.pl`, `createRepeatLandscape.pl` | none; run per genome as separate jobs |
| `summarize_rm.py`, plotting | none |
| sample-level independence | all per-sample rules run concurrently across samples |

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
per-sample chunk into the final `{outdir}/summary/*.tsv` tables and
computes `arm_concordance.tsv`. This keeps `summarize_rm.py` itself focused
on one (arm, sample) at a time (matching the spec's description of what
it does) rather than overloading it with cross-run aggregation.

## Known corrections applied vs. the spec's exact rule text (see plan for full rationale)

- `cd-hit-est` is run with `-r 1` (both strands) — the spec's command
  omitted it, which would silently under-merge reverse-complement
  duplicates between two independently-run samples (RepeatModeler2 gives
  no guarantee two samples' assemblies will report the same family on the
  same strand).
- `build_db` runs inside the same Singularity image as `repeatmodeler`,
  not the conda RepeatMasker env, to avoid a version mismatch between the
  database writer and reader.
- `library_membership.tsv` additionally reports per-cluster
  `label_agreement`/`distinct_labels`, since cd-hit-est keeps one arbitrary
  representative per cluster and silently discards the rest — if two samples' independent RepeatClassifier calls disagreed on a merged
  family's `Class/Family`, that's now visible instead of silently biasing
  `class_composition.tsv` toward whichever label cd-hit happened to keep.

## Provenance

`{outdir}/summary/provenance.txt` records RepeatMasker/RepeatModeler
versions, the FamDB release info, cd-hit version, which Dfam
partition(s) each sample's RepeatModeler container could see, the pinned
container tag, each sample's genome fingerprint, the LTR tool versions
(genometools, LTR_retriever, ltr_finder; vendored script commits are in
`workflow/vendor/README.md`), LTR candidate counts and window-timeout
skipped bp, any `curated_override` in effect, and a full `config.yaml`
snapshot.
