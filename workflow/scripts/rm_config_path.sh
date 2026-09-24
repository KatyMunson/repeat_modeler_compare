#!/bin/bash
# Print one RepeatModeler configuration value (e.g. GENOMETOOLS_DIR,
# CDHIT_DIR, LTR_RETRIEVER_DIR) as the RepeatModeler inside the current
# container resolves it. In dfam/tetools these tools are installed under
# /opt but are NOT on PATH (LTR_retriever, gt, cd-hit) -- RepeatModeler
# finds them through RepModelConfig.pm, and so do we, rather than
# hard-coding container paths.
#
# Usage: rm_config_path.sh <KEY>
set -euo pipefail
key="$1"
rm_bin="$(command -v RepeatModeler || true)"
if [ -z "$rm_bin" ]; then
    echo "[rm_config_path] RepeatModeler not on PATH -- run this inside the tetools container" >&2
    exit 1
fi
rm_dir="$(dirname "$(readlink -f "$rm_bin")")"
value="$(perl -I"$rm_dir" -MRepModelConfig -e \
    'RepModelConfig::resolveConfiguration({}); print $RepModelConfig::configuration->{$ARGV[0]}->{value} // ""' \
    "$key" 2>/dev/null || true)"
if [ -z "$value" ]; then
    echo "[rm_config_path] could not resolve $key from $rm_dir/RepModelConfig.pm" >&2
    exit 1
fi
echo "$value"
