# Spec: genome-wide NUMT calls and structural DNA-TE discovery

For a Code session working in `repeat_modeler_compare`. Two independent features; build and
merge them separately (NUMT first: small, no remasking).

Read first: `README.md` sections "Why LTR discovery runs outside RepeatModeler", "Satellite,
rDNA and mito cross-check", "Discovery summary"; `workflow/scripts/merge_families.py`
docstring. Follow repo conventions: conda envs under `workflow/envs/`, `exec > {log} 2>&1`,
`set -euo pipefail`, every config key documented in `config.yaml` and the README, report-only
by default.

## Feature A — genome-wide NUMT calls

### Why
The cross-check only labels a *library family* `Other/NUMT` when its consensus matches the
mitogenome (`satellite_consensus_blast` → `library_vs_refs.tsv`). Single-copy NUMTs never
become RepeatModeler families, so most NUMTs are not reported anywhere. The paper needs per-
locus NUMT calls, and a later base-level reconciliation will give NUMT bp priority over TE
labels.

### Inputs
- `satellite_crosscheck.external_annotations`, column `mitohifi_fa` (already parsed into
  `SATX_ANNOT`; `.`/`NA` = not provided). Samples without one are skipped, with a note.
- The prepped genome and `name_map.tsv` / `fingerprint.tsv` from `prep_genome` /
  `genome_fingerprint` (report coordinates in the original contig names, as other outputs do).

### Method (new rule `numt_calls`, per sample; new script `workflow/scripts/numt_calls.py`)
1. Mitogenome doubled (seq+seq) so hits across the circular origin aren't split; fold hit
   coordinates back modulo the mitogenome length.
2. `blastn -task dc-megablast` (`workflow/envs/blast.yaml`), mitogenome as query vs the nuclear
   genome db, `-evalue numt.evalue` (default 1e-5), `-dust no`, **`-max_target_seqs` set above
   the assembly's contig count** (default 500 would silently cap the number of contigs
   reported) and **no `-max_hsps`** (unlimited HSPs per contig = every NUMT on a chromosome).
   blastn rather than minimap2: old NUMTs are short and ~70–85 % identical, below what
   minimap2's seeds reliably find; BLAST/LAST is the usual choice in NUMT studies.
3. **Exclude mitochondrial contigs**: any contig where hits cover ≥ `numt.mito_contig_cov`
   (0.8) of the contig at ≥ `numt.mito_contig_id` (98) %. List them in the log and in the
   summary; never call NUMTs on them.
4. Drop duplicate hits from the doubled query (same genome interval).
5. Two levels of call, both written:
   - **hit-level** NUMTs (one per HSP after step 4) — the primary, unambiguous count;
   - **compound** NUMTs: adjacent hits on the same contig and strand joined only when the
     *mito* coordinates continue — the next hit starts within `numt.mito_merge_gap` (100 bp)
     of where the previous one ended (gap or overlap; origin wrap allowed) — and the *genome*
     gap is ≤ `numt.merge_gap` (default 500 bp). That is one NUMT split by an indel or a later
     insertion, not two separate insertions. Record the genome gap and what fills it
     (shared-arm RepeatMasker class, via `intervals.py`).
   - Make the default merge gap a measured choice: write a histogram of genome gaps between
     mito-colinear neighbours (`numt_gap_hist.tsv`). If gaps cluster at TE lengths (e.g. 1–6
     kb, filled by one TE), raise `merge_gap` to cover that and say so in methods. Don't start
     from 2 kb without that evidence.
   Each compound NUMT keeps: span, aligned bp, bp-weighted identity, mito start/end (and
   whether it wraps), number of hits, largest gap and its fill.
6. Keep hit-level and compound NUMTs with aligned bp ≥ `numt.min_len` (default 50).

### Outputs
| file | content |
|---|---|
| `{outdir}/numt/{sample}/numt_hits.bed` | hit-level calls: contig, start, end, id, aligned_bp, strand, identity, mito_start, mito_end |
| `{outdir}/numt/{sample}/numts.bed` | compound calls: as above + n_hits, max_gap, gap_fill |
| `{outdir}/numt/{sample}/numt_gap_hist.tsv` | genome gaps between mito-colinear neighbouring hits (for choosing `merge_gap`) |
| `{outdir}/numt/{sample}/mito_contigs.tsv` | contigs excluded as mitochondrial, with cov and identity |
| `{outdir}/summary/numt_summary.tsv` | per sample: n NUMTs, total bp, % non-N, median identity, length quantiles, bp overlapping shared-arm `.out` hits by class (`intervals.py` for the overlap) |

### Wiring
- Target rule `numt_only`; included in `all` when any sample has `mitohifi_fa`.
- Config block `numt:` with the keys above; resources block `numt_calls`.
- Report-only: nothing relabelled, nothing remasked.

### Tests / acceptance
- Unit tests for the origin-wrap folding and the merge logic (synthetic hit tables).
- Toy genome: plant 3 mitogenome fragments (one across the origin, one reverse-complemented,
  one 85 % identity) plus one full mito contig. Also split one fragment with a planted 3 kb TE
  copy. Expect 3 compound NUMTs at
  `merge_gap: 5000` (4 at the default 500), the origin-spanning one as a single call, the mito
  contig excluded.
- Real samples: the families the cross-check already calls `Other/NUMT` must overlap called
  NUMTs (sanity check; print the overlap in the log).

## Feature B — structural TIR / Helitron discovery (pilot, then optional integration)

### Why
RepeatModeler's RECON/RepeatScout rounds model TIR elements poorly and Helitrons hardly at
all. The LTR side pipeline already adds structure-based LTR families; this is the DNA-TE
analog. **Adding families changes the library, so both RepeatMasker arms and everything
downstream rerun.** Build the pilot first and only integrate if it adds real bp.

### Stage B1 — pilot (report-only, no library change)
- Rule `dna_te_candidates` per sample: run TIR and Helitron structural callers on the prepped
  genome. Preferred: EDTA's modules run standalone (TIR-Learner, HelitronScanner via
  `EDTA_raw.pl --type tir|helitron`), in their own conda env/container. Reuse
  `ltr_group_genome`'s bp-balanced groups (whole scaffolds) for SGE parallelism; check whether
  each caller needs per-window timeouts on the hagfish assemblies (the LTRharvest stall
  lesson) and log skipped regions the same way (`skipped_windows.tsv`).
- Rule `dna_te_pilot`: cluster candidates (cd-hit-est, `library.cdhit` 80/80, both strands)
  against the sample's **own** `-families.fa` and the shared library. Report in
  `{outdir}/summary/dna_te_pilot.tsv` per sample × type (TIR, Helitron): candidates, clusters,
  clusters matching an existing family (and that family's current Class/Family), unmatched
  clusters, and the bp unmatched clusters would mask (quick RepeatMasker of the genome with
  just the unmatched consensi, or blastn coverage as a cheaper proxy — pick one, say which).
- Decision rule for the user (print it in the summary): integrate if unmatched clusters add
  ≥ `dna_te.min_new_pct_non_n` (default 0.5 %, same scale as the round-novelty stopping rule)
  of the genome, or reveal a class absent from the library.

### Stage B2 — integration (only after the pilot says so; gated by `dna_te.enabled`)
- Build consensi + Stockholm seeds for unmatched clusters (RepeatClassifier needs both): MAFFT
  of up to N members per cluster + majority consensus, or RepeatModeler's `Refiner` the way
  `LTRPipeline_from_scn` does. Prefer Refiner for consistency with the LTR families.
- Family IDs: `tir-1_family-N`, `hel-1_family-N`.
- Extend `merge_families.py`: structural families (ltr, tir, hel) all have the precedence LTR
  families have now (a cluster with a structural member drops its `rnd-*` members). Document
  the change in the docstring and README.
- **Every place that assumes only `rnd-`/`ltr-` prefixes must learn the new ones**:
  `merge_families.py` (`CLSTR_MEMBER` regex), `library_membership.py` (`sources`),
  `round_saturation.py` (`LTR` regex / buckets), `round_novelty.py`, `discovery_summary.py`
  (`merged_ltr_families`; add `merged_tir_families`, `merged_hel_families`),
  `curation_candidates.py` / `family_neighbors.py` (`LTR_PIPELINE` regex — decide per script
  whether "structural" or strictly LTR is meant). Grep `ltr-` and `_ltr-` to catch the rest.
- Methods text: add to the README's "Relative to a stock run, things differ" list.
- Outputs into a new `outdir` (as `results_v2`/`results_v3` were), so the current results
  stay untouched.

### Tests / acceptance
- B1: dry run schedules only the new rules (`--rerun-triggers mtime -n`); pilot table produced
  for both hagfish samples.
- B2: `discovery_summary.tsv` per-sample family counts still reconcile with the `.clstr`
  (existing check); unit test for `merge_families.py` with tir/hel members; `round_novelty`
  and `library_membership` show the new sources.

## Out of scope here
Base-level reconciliation of satellite / TE / NUMT / rDNA labels and the combined figures
(will live in the integrate step); telomere (tidk) and tRNA/5S tracks.
