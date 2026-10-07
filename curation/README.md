# Curated elements

Full-length consensi of elements curated by hand (README "Curating a family",
"Naming curated elements"). The curated-families table that maps library
families to these elements (classify.curated_families) is built from blastn
of these consensi against the library: `family_groups.py members` per element,
then `family_groups.py resolve` across related elements.

| element | class | built from | notes |
|---|---|---|---|
| Gypsy-N1_Esto | LTR/Gypsy, non-autonomous | Esto_rnd-1_family-332 extended (18 copies) | was "EstoLTR332". ~8.1 kb, LTRs ~1.7 kb, no ORFs; GC-rich tandem repeats replace gag-pol. Shares PBS, PPT, LTR termini/U5, 5' leader and the post-pol 3' region with Gypsy-1_Esto. Inserts in (TA)n. |
| Gypsy-1_Esto | LTR/Gypsy, autonomous | Esto_rnd-5_family-3054 extended (14 copies) | ~11.8 kb, LTRs ~1.6 kb; gag (681 aa) + pol (1272 aa; PR, RT, RH, INT gmr1-like). Same PBS (tRNA primer) and PPT as Gypsy-N1_Esto. Inserts in (TA)n. |

`hagfish_curated.tsv` is the curated-families table for the hagfish run
(`classify.curated_families: curation/hagfish_curated.tsv`): blastn of both
consensi against `results_v3/classify/library_consensi.fa`, `members` per
element (`--class-family LTR/Gypsy`, `--ltr-len` 1695 / 1608), `resolve`; the
four Mlim LTR-pipeline families matching Gypsy-1_Esto (Mlim's own Gypsy
relatives, already LTR/Gypsy, < 0.03 Mb) were removed.

### Satellites, rDNA and other relabels (added to hagfish_curated.tsv)

Satellite groups are named after the satellite pipeline's harmonized motifs
(`harmonized_repeatmasker_lib.fasta`, `<SP>_SAT<unit>_<letter>`): each tandem
library family was scored against every motif tiled to the consensus length,
in 120 bp windows (>= 70% identity per window; `family_cov` = consensus
fraction covered, `pident` = median window identity).
- `part = .`: the consensus is the satellite (>= ~0.5 coverage).
- `part = partial`: a satellite segment plus other sequence (RepeatModeler
  chimeras such as Mlim_rnd-4_family-317 = MLI_SAT33_a + ~1.4 kb non-coding;
  no REXdb/GyDB domains at E <= 1, so not a Gypsy).
- `SAT-<n>_<sample>`: satellite by the cross-check (TideCluster arrays) but no
  harmonized motif matched; provisional names to reconcile with the
  satellite pipeline.
- `rDNA_45S`: ribotin rDNA matches (high).
- `TRC1-junction_Mlim` (Mlim_rnd-5_family-3676, was LTR/Gypsy): TRC_1 array
  edge + CR1 fragment + family-99 post-gag segment, no coding domains ->
  Unknown.

Rows the satellite cross-check now reproduces (harmonized motif evidence,
`satellite_crosscheck.harmonized_library`) were removed: it proposes them
at high confidence with the same name and `satellite_crosscheck.apply`
(`apply_confidence: [high]`) applies them from `summary/satellite_applied.tsv`.
These are Mlim_rnd-1_family-3 and Mlim_rnd-5_family-7820 (MLI_SAT32_a),
Mlim_rnd-5_family-3270 (MLI_SAT32_d), Mlim_rnd-1_family-587 (MLI_SAT1016_a),
Mlim_rnd-4_family-274 (MYX2_SAT49_a), Mlim_rnd-5_family-7244
(MLI_SAT1980_a), Esto_rnd-1_family-686 (EST_SAT179_b) and
Esto_rnd-5_family-6663 (EST_SAT41_c). Mlim_rnd-3_family-37 was also removed:
the cross-check names it MYX2_SAT25_a (25 bp motif = its 25 bp KITE founder),
replacing the earlier windowed MYX2_SAT36_a. Rows kept here are medium or
partial calls (not auto-applied), provisional SAT-n names, rDNA_45S and the
junction; if `apply` is turned off, restore the removed rows from git history.
