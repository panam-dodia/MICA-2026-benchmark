#!/bin/bash
# Run ModelAngelo (sequence-guided: model_angelo build -v MAP -pf FASTA -o OUT) on every entry in runs/ready_entries.tsv.
# Resumable: entries with a DONE file are skipped. Optional: ONLY=<PDB ID> runs one entry.
#
# Runs alongside the MICA batch, so before each entry it waits until the GPU has >= GPU_MIN_MIB free and the
# system has >= RAM_MIN_GIB available. The batches overlap, so runtimes are approximate (see notes).
# Inputs per entry: map + RCSB FASTA only (no ground truth, no AF3). Final model: <entry>/output/output.cif
#
# Usage:  bash scripts/run_modelangelo_batch.sh [entries.tsv]
set -u

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
RUNS=$RUNS_ROOT/modelangelo
LIST=${1:-$BENCH/runs/ready_entries.tsv}
TIMES=$RUNS/timings.csv
ONLY=${ONLY:-}
GPU_MIN_MIB=${GPU_MIN_MIB:-12288}
RAM_MIN_GIB=${RAM_MIN_GIB:-12}

source "$CONDA_SH"
conda activate model_angelo

mkdir -p "$RUNS"
[ -f "$TIMES" ] || echo "emdb,pdb,step,start,end,seconds,exit_code" > "$TIMES"

gpu_free() { nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1; }
ram_avail() { awk '/MemAvailable/ {print int($2 / 1024 / 1024)}' /proc/meminfo; }

wait_for_resources() {
    local waited=0
    while true; do
        local g r
        g=$(gpu_free); r=$(ram_avail)
        if [ "$g" -ge "$GPU_MIN_MIB" ] && [ "$r" -ge "$RAM_MIN_GIB" ]; then
            [ "$waited" -gt 0 ] && echo "    resources free after ${waited} s wait (GPU free ${g} MiB, RAM avail ${r} GiB)"
            return
        fi
        [ $((waited % 600)) -eq 0 ] && echo "    waiting: GPU free ${g} MiB (need ${GPU_MIN_MIB}), RAM avail ${r} GiB (need ${RAM_MIN_GIB})"
        sleep 30; waited=$((waited + 30))
    done
}

tail -n +2 "$LIST" | while IFS=$'\t' read -r N PDB CONTOUR RES RESIDUES JOBS; do
    [ -n "$N" ] || continue
    if [ -n "$ONLY" ] && [ "$PDB" != "$ONLY" ]; then continue; fi
    D=$RUNS/$N
    if [ -f "$D/DONE" ]; then echo "skip $PDB (done)"; continue; fi

    echo "=== $(date '+%F %T')  $PDB / EMD-$N  ($RESIDUES deposited residues)"
    rm -rf "${D:?}"
    mkdir -p "$D"
    cp "$BENCH/dataset/entries/$N/emd_$N.map" "$D/"
    cp "$BENCH/dataset/entries/$N/${PDB,,}.fasta" "$D/"

    wait_for_resources
    s=$(date +%s)
    model_angelo build -v "$D/emd_$N.map" -pf "$D/${PDB,,}.fasta" -o "$D/output" --device 0 < /dev/null > "$D/modelangelo.log" 2>&1
    rc=$?
    e=$(date +%s)
    echo "$N,$PDB,build,$(date -d @"$s" +%FT%T),$(date -d @"$e" +%FT%T),$((e - s)),$rc" >> "$TIMES"
    if [ $rc -eq 0 ] && [ -s "$D/output/output.cif" ]; then
        echo "    build: $((e - s)) s, model: $D/output/output.cif"
        touch "$D/DONE"
    else
        echo "    build FAILED (exit $rc) after $((e - s)) s; see $D/output/model_angelo.log and $D/modelangelo.log"
        echo "build exit $rc" > "$D/FAILED"
    fi
done
echo "=== batch finished $(date '+%F %T')"
