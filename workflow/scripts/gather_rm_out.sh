#!/bin/bash
# Adapted from KatyMunson/compare_assemblies_satellites
# (common/scripts/gather_rm_out.sh, commit 27ce4ee) for repeat_compare's
# repeatmasker scatter/gather. One real bug fixed (see the `|| true` below):
# the original's `set +o pipefail` only protects the header-seed line
# (its pipe ends in `head`, immune to grep's exit code); the for-loop's
# pipe ends in `grep -v "^There"` itself, whose "no lines matched" exit
# code (1) still trips `set -e` when a chunk has genuinely zero hits
# (exactly what a touch()-only chunk .out produces) -- reproduced locally
# against a fixture, confirmed this aborts the whole script mid-merge.
# Never surfaced upstream since their real-genome chunks apparently never
# came up fully empty; a heavily subsampled wiring-test genome split into
# many small chunks hits this routinely. Otherwise unmodified.
#
# gather_rm_out.sh -- merge N per-chunk RepeatMasker .out files into one.
#
# Shared logic, previously duplicated byte-for-byte between
# compare_meadowlark's gather_repeatmasker rule and trf_consensus's
# satellite_screen (now stage 03 motif_evaluation) gather_repeatmasker_genome_lib
# rule. RepeatMasker .out has a fixed 3-line header block, plus an optional
# "There were no repetitive sequences detected" sentinel line for a chunk
# with zero hits. Not a plain concatenation -- that would repeat the header
# block once per chunk: keep the first input's first 3 lines (its header,
# sentinel-stripped) as the seed, then append every input's data rows (line
# 4 onward, sentinel-stripped).
#
# Usage: gather_rm_out.sh <output_path> <chunk1.out> [chunk2.out ...]
set -euo pipefail

out="$1"
shift

set +o pipefail  # grep -v legitimately returns nonzero when a chunk has nothing left after filtering
first="$1"
grep -v "^There" "$first" | head -n 3 > "$out"
for f in "$@"; do
    tail -n +4 "$f" | grep -v "^There" || true
done >> "$out"
set -o pipefail
