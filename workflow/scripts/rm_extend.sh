#!/usr/bin/env bash
# Run, recover or extend RepeatModeler (rounds only) for one sample, then copy
# the cumulative consensi.fa / families.stk out. Called by rule repeatmodeler
# from the sample's repeatmodeler/ directory (where BuildDatabase wrote the db).
#
# usage: rm_extend.sh SAMPLE THREADS EXTRA_ROUNDS "EXTRA_ARGS" FINGERPRINT FINGERPRINT_PY \
#                     OUT_CONSENSI OUT_STK OUT_PROV
#
# Target: T = 5 + EXTRA_ROUNDS rounds (RepeatModeler 2.0.9: round 1
# RepeatScout, rounds 2-5 RECON at 10/30/90/270 Mb; every extra round
# samples another 270 Mb).
#
# How RepeatModeler 2.0.9 recovery works (RepeatModeler source, lines
# 678-790 and 862-880 / 1444-1475):
#  - -recoverDir takes the highest round-N/ with a non-empty consensi.fa as
#    the last good round h, and reruns round h+1 -- but only if a round-(h+1)/
#    directory exists (empty is fine); otherwise it prints "appears to contain
#    a successful run" and exits.
#  - New families are appended to the top-level consensi.fa / families.stk,
#    which must therefore stay in place.
#  - When resuming at a round already at the 270 Mb cap, -numAddlRounds N
#    runs N+1 rounds (the resumed round is not counted). To end at round T
#    from h >= 5, pass T-h-1. A resume below the cap (h+1 <= 5) behaves like
#    a fresh run: pass EXTRA_ROUNDS.
#
# Extending a finished run never touches it: RM_x is copied to RM_x.ext and
# the copy is extended; rounds 1-5 of the copy must stay byte-identical.
set -euo pipefail

SAMPLE=$1 THREADS=$2 EXTRA=$3 EXTRA_ARGS=$4 FP=$5 FP_PY=$6 OUT_CONS=$7 OUT_STK=$8 OUT_PROV=$9
T=$((5 + EXTRA))

highest() {  # highest round-N with a non-empty consensi.fa
    local d=$1 h=0 i=1
    while [ -s "$d/round-$i/consensi.fa" ]; do h=$i; i=$((i + 1)); done
    echo "$h"
}
addl_for() {  # -numAddlRounds to reach T when resuming after round $1
    local h=$1
    if [ $((h + 1)) -le 5 ]; then echo "$EXTRA"; else echo $((T - h - 1)); fi
}
same_rounds() {  # rounds 1-5 of $1 and $2 byte-identical
    for i in 1 2 3 4 5; do
        if ! cmp -s "$1/round-$i/consensi.fa" "$2/round-$i/consensi.fa"; then
            echo "[ERROR] round-$i/consensi.fa differs between $1 and $2; refusing to use $2."
            echo "[ERROR] Delete $2 to start the extension again from $1."
            exit 1
        fi
    done
}
recover() {  # recover/extend directory $1 up to round T
    local d=$1 h addl
    h=$(highest "$d")
    if [ "$h" -ge "$T" ]; then return 0; fi
    addl=$(addl_for "$h")
    mkdir -p "$d/round-$((h + 1))"
    echo "[INFO] $d: rounds 1-$h done; running round(s) $((h + 1))-$T (-numAddlRounds $addl)"
    # shellcheck disable=SC2086
    RepeatModeler -database "$SAMPLE" -threads "$THREADS" -recoverDir "$d" -numAddlRounds "$addl" $EXTRA_ARGS
    h=$(highest "$d")
    if [ "$h" -lt "$T" ]; then
        echo "[ERROR] $d ends at round $h after recovery; expected $T."
        exit 1
    fi
}

RM_DIR=$(ls -d RM_* 2>/dev/null | grep -v '\.ext$' | head -n1 || true)
USE=""
if [ -z "$RM_DIR" ]; then
    cp "$FP" rm_run.fingerprint.tsv
    ADDL=""
    if [ "$EXTRA" -gt 0 ]; then ADDL="-numAddlRounds $EXTRA"; fi
    # shellcheck disable=SC2086
    RepeatModeler -database "$SAMPLE" -threads "$THREADS" $ADDL $EXTRA_ARGS
    RM_DIR=$(ls -d RM_* | grep -v '\.ext$' | head -n1)
    USE=$RM_DIR
else
    if [ ! -s rm_run.fingerprint.tsv ]; then
        echo "[ERROR] $RM_DIR exists but rm_run.fingerprint.tsv doesn't -- it wasn't started by this rule."
        echo "[ERROR] Move $RM_DIR out of $(pwd) (or delete it) and rerun."
        exit 1
    fi
    python3 "$FP_PY" check --expected rm_run.fingerprint.tsv --observed "$FP" || {
        echo "[ERROR] $RM_DIR was started on a different genome than the current prepped FASTA."
        echo "[ERROR] Refusing to -recoverDir it. Move $RM_DIR and rm_run.fingerprint.tsv away and rerun."
        exit 1
    }
    if [ ! -s "$RM_DIR/consensi.fa.classified" ]; then
        # an interrupted run: finish it in place (up to T rounds)
        echo "[INFO] Recovering interrupted RepeatModeler run $RM_DIR"
        recover "$RM_DIR"
        USE=$RM_DIR
    elif [ "$(highest "$RM_DIR")" -ge "$T" ]; then
        echo "[INFO] $RM_DIR already has $(highest "$RM_DIR") rounds (target $T); reusing it"
        USE=$RM_DIR
    else
        EXT="$RM_DIR.ext"
        if [ ! -d "$EXT" ]; then
            echo "[INFO] copying finished $RM_DIR to $EXT to extend it (the original is not modified)"
            cp -a "$RM_DIR" "$EXT"
        fi
        same_rounds "$RM_DIR" "$EXT"
        recover "$EXT"
        same_rounds "$RM_DIR" "$EXT"
        USE=$EXT
    fi
fi

test -s "$USE/consensi.fa" && test -s "$USE/families.stk"
cp "$USE/consensi.fa" "$OUT_CONS"
cp "$USE/families.stk" "$OUT_STK"
{
    echo "=== RepeatModeler version ==="
    RepeatModeler -version 2>&1 || echo "RepeatModeler -version failed"
    echo "=== run: rounds only (no -LTRStruct), extra_rounds=$EXTRA (target $T rounds; $(highest "$USE") present), extra_args='$EXTRA_ARGS', dir $USE ==="
    if [ "$USE" != "$RM_DIR" ]; then echo "extended copy of $RM_DIR (rounds 1-5 verified identical)"; fi
    echo
    echo "=== famdb.py info (as seen by RepeatClassifier inside this container) ==="
    famdb.py info 2>&1 || echo "famdb.py info failed or not found in container"
} > "$OUT_PROV"
