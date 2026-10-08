# Manual curation: state and to-do list (hagfish, Esto / Mlim)

Working notes for the next curation pass. Family IDs are library IDs
(`<sample>_rnd-N_family-M` from RepeatModeler, `<sample>_ltr-1_family-M`
from the LTR pipeline). Numbers refer to the shared arm
(`results_v3/summary_shared_only/`) as of the last `family_neighbors` run
(commit `cd61760`: orientation-aware, binomial test on pooled LTR partners).

## Where things stand

- `curation/hagfish_curated.tsv` (111 rows) holds:
  - Gypsy-N1_Esto / Gypsy-1_Esto members;
  - satellites named after harmonized motifs, plus provisional `SAT-n` names;
  - rDNA_45S;
  - the TRC_1 junction;
  - 35 neighbour-placed families (safe set from the first neighbour run).
- `summary/satellite_applied.tsv` applies the high-confidence satellite
  cross-check calls (9 former curated rows now come from there).
- Last `family_neighbors` run:
  - 2,467 Unknown / `<Order>/Unknown` families tested;
  - 521 proposals (9 high), covering ~0.98 M of 4.9 M Unknown copies (20%);
  - 157 Unknown chains.
- Unknown still ~566 Mb in Esto (22.9%) and ~601 Mb in Mlim (15.7%)
  before applying any neighbour calls.

## 1. Next additions to the curated table (from the last neighbour run)

Take from `summary_shared_only/family_neighbors_proposals.tsv`:
- the 9 high-confidence calls;
- the 79 `internal_of_LTRs` calls (LTR pieces on both sides);
- Mlim_rnd-5_family-7567 → LTR/Gypsy, the 3' piece of Mlim_rnd-1_family-50
  (42% of its copies, opposite strand).

Keep as evidence only, not labels, for now:
- the 245 `next_to_LTRs` calls;
- the one-sided `5prime_of` / `3prime_of` medium calls;
- all partner-anchored and propagated calls.

Check before accepting:
- Esto_rnd-1_family-361 → LTR/ERVK (anchored by Esto_rnd-4_family-1923, 52%).
  Earlier runs also tied it to a hAT, and it carries a tandem note.
- Esto_rnd-1_family-29 → DNA/Maverick via Esto_ltr-1_family-338. An LTR-pipeline
  family labelled Maverick is suspect (see section 4).
- Esto_rnd-2_family-60 (17,300 copies): partner Esto_rnd-1_family-284 is
  labelled `RNA`, so the label was not passed on. Check what 284 is.

## 2. The large Mlim Gypsy element (rnd-1_50 group)

- Members so far:
  - Mlim_rnd-1_family-50 (LTR/Gypsy);
  - Mlim_rnd-5_family-7567;
  - Mlim_ltr-1_family-873 (high, reverse-strand 3' piece of 50);
  - Mlim_rnd-1_family-8 (lost its call under the stricter test, but 50 is its
    top 5' partner).
- Curate like Gypsy-N1_Esto:
  1. `family_profile.py` on 50 and 7567;
  2. extend to full-length copies;
  3. find the LTRs, PBS and PPT;
  4. `family_groups.py members` against the assembled element;
  5. name it (`Gypsy-<n>_Mlim`).

## 3. Mlim composite block (8 / 544 / 353 / 221 / 70 / 772 / 50)

- These families sit next to each other in fixed orientations across
  thousands of copies:
  - Mlim_rnd-1_family-8, -544, -353, -221, -70;
  - the DNA/Academ-1 family Mlim_rnd-5_family-772;
  - the Gypsy Mlim_rnd-1_family-50.

  That includes a DNA transposon, so it is not one LTR element. Most likely a
  multi-kb composite unit amplified as a block. That fits repeat-rich
  germline-restricted chromatin, but this is not yet tested.
- Their pooled LTR calls did not survive the binomial test. Do not label them
  `LTR/Unknown`.
- To do:
  1. Write `block_loci.py`: find loci where several of these families co-occur
     in the `.out`, and export coordinates plus ±25 kb of sequence.
  2. Dot-plot a handful of loci to get the unit length and structure.
  3. Check whether the block is Mlim-specific and, given the germline/somatic
     question, whether it is enriched on the eliminated fraction (needs the
     depth data, see section 9).
- Related large Unknowns, also candidates for block membership:
  - Mlim_rnd-1_family-26 (37,700 copies; chain-1 with Esto_rnd-1_family-319;
    matched satellite motif MLI_SAT955_a at 97.8% earlier, tandem_frac 0.15);
  - Mlim_rnd-1_family-462;
  - Mlim_rnd-1_family-407 (chain-2);
  - Mlim_rnd-3_family-628 (chain-4).

## 4. Other families flagged during this work

- **Mlim_rnd-5_family-769:** 20 kb consensus, ~140 Mb, labelled LINE by the
  reclassification.
  - 20 kb is too long for a LINE: check for a chimera or Maverick.
  - Mlim_rnd-4_family-104 (20,500 copies) is its 3' piece in 60% of copies.
  - Mlim_rnd-1_family-229 joins through 104.
  - Run `family_profile.py` first.
- **Maverick wave (Esto):** several families labelled DNA/Maverick, some from
  the LTR pipeline (e.g. Esto_ltr-1_family-338). Recheck with the
  integrase-like / Maverick rule in `reclassify_unknown.py` in mind.
- **Esto_rnd-3_family-410:** labelled DNA/PiggyBac, ~35 Mb, matched satellite
  motif EST_SAT21_a (0.66 at 98%). Check whether it is a TE-derived satellite
  or a real PiggyBac.
- **Esto_rnd-5_family-232:** LINE/RTE-BovB label, but a Gypsy domain hit,
  16 Mb. Chimera?
- **Young Unknowns (Esto):** ~156 Mb of Esto Unknown is < 4% divergent.
  - Largest: Esto_rnd-1_family-361, -29, -389, -46, -486; Esto_rnd-2_family-361.
  - Near-identical copies make these the easiest to curate (extend, find ends,
    TSDs).
  - Chain-7 (Esto_rnd-1_family-389, -36, -758) is likely one element.
- **TE-derived "satellite" motifs** to report back to the satellite pipeline.
  These are TE consensi matching harmonized motifs at ≥ 94% over most of
  their length (satellite cross-check, medium):

  | family | label | motif |
  |---|---|---|
  | Mlim_rnd-1_family-99 | L1-Tx1 | MLI_SAT1365_a |
  | Mlim_rnd-5_family-816 | L2 | MLI_SAT2379_a |
  | Mlim_rnd-5_family-97 | SINE/tRNA | MLI_SAT29_i |
  | Mlim_rnd-4_family-2271 | RTE-BovB | MLI_SAT443_a |
  | Mlim_ltr-1_family-158 | — | MLI_SAT71_g |
  | Mlim_ltr-1_family-362 | — | MLI_SAT31_e |
  | Esto_rnd-5_family-5412 | PIF-Harbinger | MYX2_SAT17_c |
  | Mlim_rnd-5_family-1180 | TcMar-Fot1 | MLI_SAT41_g |

## 5. Satellite questions for the satellite pipeline

- Esto_rnd-5_family-1124 (SAT-1_Esto):
  - likely a diverged EST_SAT225_a;
  - same 225 bp unit and same TRC_3 as Esto_ltr-1_family-61;
  - but only 13% of the consensus matches the motif at ≥ 75%.
- Mlim_rnd-3_family-37 is now MYX2_SAT25_a (25 bp KITE founder = 25 bp motif),
  replacing the earlier windowed MYX2_SAT36_a. Confirm.
- Provisional names SAT-1..7_Mlim and SAT-2_Esto: no harmonized motif matched
  well enough. Ask for IDs.
- Dfam satellites called by the cross-check (SAT-9_DR, MSAT-4_DR, HSATII, ...):
  harmonized IDs?
- Misses below threshold, already in the curated table:
  - Mlim_rnd-3_family-866 (MLI_SAT187_a, 0.47 coverage);
  - Esto_rnd-5_family-3236 (EST_SAT51_c).

## 6. TEtrimmer (target `tetrimmer_only`)

- Not run yet. `tetrimmer_select` already skips:
  - families resolved by curated / applied labels;
  - families with high-confidence neighbour calls;
  - tandem families.
- Natural first input: the chain hubs and the large uncalled Unknowns.
  - Chain hubs: Mlim_rnd-1_family-26, -407, Mlim_rnd-3_family-628,
    Esto_rnd-1_family-389, Mlim_rnd-1_family-544, -353.
  - Large uncalled: Esto_rnd-1_family-457 (32,200 copies), Mlim_rnd-1_family-48
    (78,500 copies).
- Before running:
  - pin the conda env;
  - set `tetrimmer.pfam_dir` if compute nodes have no internet.
- Possible improvement: pass chain hubs only (TEtrimmer's extension reaches
  across neighbouring pieces), not every member.

## 7. Library release (deferred: `export_library`)

- For publication, release an oriented, renamed copy of the library. Keep the
  masking library unchanged so `.out` coordinates stay valid.
- Orientation cues, strongest first:
  1. coding strand (TEsorter / ORFs);
  2. LTR structure (PBS / PPT);
  3. 5'-truncation coverage profile + poly(A) tail;
  4. propagation from oriented partners (same/opposite strand from
     `family_neighbors`);
  5. harmonized motif orientation for satellites.

  No cue → leave as is and flag `unoriented`.
- Rename to curated names (scheme A), and publish curated elements as their
  assembled consensi (`curation/elements/`).
- Write `library_release.tsv` with columns: old ID, new name, class,
  orientation, cue, reverse-complemented yes/no.

## 8. NUMT follow-ups (Mlim; `numt_only`, `numt_gap_check.py`)

Context: `results_v3/numt/Mlim/numts.bed` after the BLAST-name fix (the
assembly's own mitogenome NC_002639.1 excluded as the mito contig). The 20
compound calls with >= 5 kb aligned were checked for flanking assembly gaps
(`workflow/scripts/numt_gap_check.py`, output
`results_v3/numt/Mlim/numt_gap_check.tsv`):
- none has a gap or contig end within 5 kb on both sides (chance 0.0001),
  and 17 / 20 have no gap within 10 kb on either side, so they are not
  scaffolded mito contigs. Keep them as NUMTs;
- 2 / 20 have a gap on one side within 5 kb (10 % vs 1.5 % by chance;
  small numbers, and repeat-rich regions are harder to assemble).

To check:
- **`Mlim_numt_485`** (NW_027149408.1:470775-478247, 7.5 kb, 97.0 %, minus
  strand) ends exactly at the scaffold end. The scaffold may end inside a
  NUMT, or the assembler may have attached mito sequence. Dot plot of the
  last ~20 kb of NW_027149408.1 against the mitogenome. The left side is
  continuous for 470 kb, so it is not a mito contig on its own.
- **One-sided 500-N gaps** (NCBI's standard gap size, so assembly joins):
  - `Mlim_numt_120` (NC_090426.1:127037945-127044612): gap 2.6 kb to the left;
  - `Mlim_numt_388` (NW_027146820.1:33081-38642, a 70 kb scaffold): gap
    1.7 kb to the right.
- **Two identity groups** among the long calls:
  - about 96-97 % (9 calls): mostly single NUMTs on chromosomes
    (NC_090422 / 425 / 426 / 432). These look like independent recent
    insertions;
  - about 92-93 % (11 calls), clustered:
    - NW_027145999.1: 4 calls between 98 and 297 kb;
    - NW_027147682.1: 3 calls within 30 kb (`numt_454` and `numt_455`
      overlap by 1.4 kb on the same strand but were not merged, so their
      mito coordinates don't continue);
    - NW_027145878.1: 2 calls 20 kb apart;
    - NC_090430.1: 2 calls 30 kb apart near the chromosome start.
  - Hypothesis: one old insertion, later duplicated (segmental duplication
    of a NUMT region), rather than 11 events. Test: all-vs-all blastn of
    the 92-93 % copies. If they match each other at >= 98 % but the
    mitogenome only at 92 %, they were copied after insertion. This matters
    for counting insertion events, not for NUMT bp.
- **`merge_gap`:** still to choose from `numt_gap_hist.tsv` after the rerun
  with real gap fills (Mlim has ~13 mito-colinear pairs with 4-8 kb genome
  gaps; Esto has a single pair, > 50 kb apart).

## 9. Other open threads

- `family_profile.py` multi-family mode (several families profiled as one
  element), promised earlier.
- Germline elimination / depth analysis: waiting on HiFi BAM / mosdepth output
  and any somatic data.
  - Needed to test whether the Mlim composite block and the large Mlim
    satellites sit in eliminated DNA.
  - This is a confound for every Esto vs Mlim comparison: Mlim testis is mostly
    germline, Esto ovary roughly 50/50.
- Reading:
  - Havecker et al. 2004 (*Genome Biology* 5:225);
  - Goubert et al. 2022 (*Mobile DNA* 13:7);
  - Wicker et al. 2007 (*Nat Rev Genet* 8:973).
