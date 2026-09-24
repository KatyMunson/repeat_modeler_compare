# Vendored third-party code

Everything here is a copy of an upstream file, plus a small, commented
local patch. The pristine copies went in first, in their own commit
("Vendor upstream ... (unmodified)"), so `git diff <that commit> -- workflow/vendor`
shows exactly what repeat_compare changed.

| Path | Upstream | Version | License |
|---|---|---|---|
| `LTR_HARVEST_parallel/` | https://github.com/oushujun/LTR_HARVEST_parallel | commit `c3c9b3c` (v1.3) | MIT |
| `LTR_FINDER_parallel/` | https://github.com/oushujun/LTR_FINDER_parallel | commit `f1036ca` (v1.4) | MIT |
| `RepeatModeler/LTRPipeline_from_scn` | https://github.com/Dfam-consortium/RepeatModeler `LTRPipeline` | tag `2.0.9` | OSL-2.1 (`RepeatModeler/LICENSE`) |

The `ltr_finder` binary and its tRNA database, which LTR_FINDER_parallel's
repo bundles, are **not** vendored. The `ltr_finder_group` rule installs
`ltr_finder` from bioconda through `workflow/envs/ltr_finder.yaml` instead.

## Local patches

### LTR_HARVEST_parallel and LTR_FINDER_parallel (same three changes each)
1. **`-harvest_args` / `-finder_args`**: overrides the hard-coded tool
   parameters. repeat_compare passes RepeatModeler's own ltrharvest
   parameters, which are the `gt ltrharvest` defaults
   (`ltr_discovery.ltrharvest_args: ""`).
2. **`-timeout_log FILE`**: appends `piece<TAB>bp<TAB>salvaged|skipped` for
   every piece killed by the `-time` limit. `ltr_gather` converts piece
   names back to genome intervals and reports the union as skipped bp.
3. **The salvage re-run actually splits.** Upstream re-ran a timed-out piece
   with `-threads 1`, which takes the tool's single-threaded branch. That
   branch ignores `-size` and `-time` and re-runs the whole piece with **no
   timeout**, so a piece that stalled once could stall forever. That's the
   failure this pipeline exists to avoid. The patch uses `-threads 2`, so
   salvage splits the piece into 50 kb windows with a 20 kb overlap, each
   with its own short timeout.

   Upstream LTR_HARVEST_parallel also never passed an overlap here. The
   default 100 kb overlap is larger than the 50 kb piece, so `cut.pl` refuses
   it. Upstream never reached that code path, because of the `-threads 1`
   issue above.

Both scripts decide "timed out" from a non-zero exit status, so they treat
any failure as a timeout. Always pass an **absolute** `-gt` / `-finder`
path: both scripts `chdir` into a work directory, and a relative tool path
then fails for every piece. That shows up as every piece "salvaged", then
"skipped".

Tested locally with genometools 1.6.5 on a synthetic genome with 12
inserted LTR elements:
- the normal run finds all 12 at genome-scale coordinates;
- a slow-`gt` shim forcing every 100 kb piece to time out yields the
  identical 12 through salvage;
- a shim that also stalls the 50 kb salvage pieces logs them all as
  `skipped`.

### RepeatModeler/LTRPipeline_from_scn
RepeatModeler 2.0.9's `LTRPipeline` with two changes, both commented in
the file:
1. `-inscn <file>` replaces the internal whole-genome `gt suffixerator` +
   `gt ltrharvest` step (the single process that stalled). The supplied
   `.scn` gets the same line filtering `runLtrHarvest` applies to real
   ltrharvest output. LTR_retriever (`seqN` renaming, `-noanno`), MAFFT,
   NINJA and Refiner then run unchanged.
2. RepeatModeler's Perl modules and `Refiner` are found through `$RM_DIR`,
   or the directory of the `RepeatModeler` on `PATH`, instead of `FindBin`.
   This copy lives outside the RepeatModeler install.

The `.scn`'s `seq-nr` column must be the 0-based index of each sequence in
the FASTA passed to it, exactly as a whole-genome `gt ltrharvest` would
report it. `workflow/scripts/normalize_scn.py` guarantees that, and checks
every coordinate against the sequence length.
