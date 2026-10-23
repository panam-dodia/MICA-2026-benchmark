#!/bin/bash
# Run MICA on every entry in runs/ready_entries.tsv, one at a time. Resumable: entries with a DONE file are skipped.
#
# Per entry this runs the same three commands, with the same arguments, as MICA's MICA_pipeline.sh v1.2:
#   1. utils/process_AF3_results.py   2. utils/dock_in_map.py (Phenix)   3. run.py --run_pulchra --run_phenix
# The only addition is "-o <entry>/output" so each entry keeps its own output folder.
# run.py starts phenix.real_space_refine in the background (modules/phenix.sh ends with '&'), so we wait for it
# before the next entry. That keeps timings honest and stops refinements from piling up.
#
# Inputs are copied from $BENCH into ~/bench/runs/mica/<EMDB>/ with no ground-truth files.
# The path must not contain "results" outside "AF3_results": process_AF3_results.py does path.replace("results", ...).
#
# Usage:  bash scripts/run_mica_batch.sh [entries.tsv]
set -u

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
MICA=$METHODS_DIR/MICA
RUNS=$RUNS_ROOT/mica
PULCHRA=$MICA/modules/pulchra304/src/pulchra
LIST=${1:-$BENCH/runs/ready_entries.tsv}
TIMES=$RUNS/timings.csv

# Needed to avoid the multiprocessing deadlock seen in the sample run (docs/PROTOCOL.md)
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
# Merizo imports torch (libgomp) before numpy; Intel MKL then refuses its default INTEL threading layer
export MKL_THREADING_LAYER=GNU
# Optional: ONLY=<PDB ID> runs just that entry (used for the single-entry test)
ONLY=${ONLY:-}
source "$CONDA_SH"
conda activate MICA

mkdir -p "$RUNS"
[ -f "$TIMES" ] || echo "emdb,pdb,step,start,end,seconds,exit_code" > "$TIMES"
cd "$MICA" || exit 1

# step <name> <command...>: run with stdin from /dev/null (so it cannot eat the entry list), log output, record time
step() {
    local name=$1; shift
    local s e rc
    s=$(date +%s)
    "$@" < /dev/null >> "$LOG" 2>&1
    rc=$?
    e=$(date +%s)
    echo "$N,$PDB,$name,$(date -d @"$s" +%FT%T),$(date -d @"$e" +%FT%T),$((e - s)),$rc" >> "$TIMES"
    echo "    $name: $((e - s)) s, exit $rc"
    return $rc
}

tail -n +2 "$LIST" | while IFS=$'\t' read -r N PDB CONTOUR RES RESIDUES JOBS; do
    [ -n "$N" ] || continue
    if [ -n "$ONLY" ] && [ "$PDB" != "$ONLY" ]; then continue; fi
    D=$RUNS/$N
    OUT=$D/output
    LOG=$D/mica.log
    if [ -f "$D/DONE" ]; then echo "skip $PDB (done)"; continue; fi

    echo "=== $(date '+%F %T')  $PDB / EMD-$N  ($RESIDUES residues)"
    OOM_BEFORE=$(awk '/^oom_kill / {print $2}' /proc/vmstat)
    rm -rf "${D:?}"
    mkdir -p "$D/AF3_results" "$OUT"
    cp "$BENCH/dataset/entries/$N/emd_$N.map" "$D/"
    cp "$BENCH/dataset/entries/$N/${PDB,,}.fasta" "$D/"
    for J in ${JOBS//,/ }; do
        mkdir -p "$D/AF3_results/$J"
        cp "$BENCH/af3/results/fold_$J/fold_${J}_model_0.cif" "$D/AF3_results/$J/"
    done
    FASTA=$D/${PDB,,}.fasta
    MAP=$D/emd_$N.map
    AF3=$D/AF3_results

    step process_AF3 python utils/process_AF3_results.py -f "$FASTA" -a "$AF3" \
        || { echo "process_AF3" > "$D/FAILED"; continue; }
    # process_AF3_results.py exits 0 even when Merizo fails, so check that it actually produced domains
    NDOM=$(ls "$D/AF3_domains"/*.pdb 2>/dev/null | wc -l)
    echo "    AF3 domains: $NDOM"
    if [ "$NDOM" -eq 0 ]; then echo "process_AF3: no domains (see mica.log)" > "$D/FAILED"; continue; fi
    step dock_in_map python utils/dock_in_map.py -m "$MAP" -c "$CONTOUR" -r "$RES" -f "$FASTA" -a "$AF3" --phenix_act "$PHENIX_ACT" \
        || { echo "dock_in_map" > "$D/FAILED"; continue; }
    step run_py python run.py -m "$MAP" -f "$FASTA" -a "$AF3" --run_pulchra --pulchra_path="$PULCHRA" \
        --run_phenix --phenix_act="$PHENIX_ACT" --resolution="$RES" --device=cuda -o "$OUT" \
        || { echo "run_py" > "$D/FAILED"; continue; }

    # wait for the detached phenix.real_space_refine (its command line contains the output path)
    s=$(date +%s); sleep 20
    while pgrep -f "real_space_refine.*$OUT" > /dev/null; do sleep 30; done
    e=$(date +%s)
    REFINED=$(ls "$OUT"/*real_space_refined*.pdb 2>/dev/null | head -1)
    echo "$N,$PDB,phenix_refine_wait,$(date -d @"$s" +%FT%T),$(date -d @"$e" +%FT%T),$((e - s)),$([ -n "$REFINED" ] && echo 0 || echo 1)" >> "$TIMES"
    echo "    phenix_refine: $((e - s)) s, refined model: ${REFINED:-MISSING}"

    # Phenix exits 0 even if its worker processes were OOM-killed, so compare the kernel OOM-kill counter
    OOM_AFTER=$(awk '/^oom_kill / {print $2}' /proc/vmstat)
    if [ "$OOM_AFTER" -gt "$OOM_BEFORE" ]; then
        echo "    WARNING: $((OOM_AFTER - OOM_BEFORE)) process(es) OOM-killed during this entry - results not trusted"
        echo "oom_kills $((OOM_AFTER - OOM_BEFORE))" > "$D/FAILED"
        continue
    fi
    touch "$D/DONE"
done
echo "=== batch finished $(date '+%F %T')"
