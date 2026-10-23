#!/bin/bash
# Run CryoAtom (CryoAtom2 v2.1.1, sequence-guided: cryoatom build --map-path MAP --protein-sequence-path FASTA --output-dir OUT) on every entry in runs/ready_entries.tsv.
# Resumable: entries with a DONE file are skipped. Optional: ONLY=<PDB ID> runs one entry.
#
# Before each entry it waits until the GPU has >= GPU_MIN_MIB free and >= RAM_MIN_GIB RAM is available (safety
# check; methods are meant to run one after another). FASTA: RCSB per-entity (one record per unique sequence, as the
# CryoAtom README asks). Language models (ESM2, RNA-FM) are pre-cached in ~/.cache/torch/hub/checkpoints.
# Inputs per entry: map + RCSB FASTA only (no ground truth, no AF3). Final model: <entry>/output/output.cif
#
# Usage:  bash scripts/run_cryoatom_batch.sh [entries.tsv]
set -u

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
RUNS=$RUNS_ROOT/cryoatom
LIST=${1:-$BENCH/runs/ready_entries.tsv}
TIMES=$RUNS/timings.csv
ONLY=${ONLY:-}
GPU_MIN_MIB=${GPU_MIN_MIB:-15360}
RAM_MIN_GIB=${RAM_MIN_GIB:-12}

source "$CONDA_SH"
conda activate CryoAtom2

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
    cryoatom build --map-path "$D/emd_$N.map" --protein-sequence-path "$D/${PDB,,}.fasta" --output-dir "$D/output" --device 0 < /dev/null > "$D/cryoatom.log" 2>&1
    rc=$?
    e=$(date +%s)
    echo "$N,$PDB,build,$(date -d @"$s" +%FT%T),$(date -d @"$e" +%FT%T),$((e - s)),$rc" >> "$TIMES"
    if [ $rc -eq 0 ] && [ -s "$D/output/output.cif" ]; then
        echo "    build: $((e - s)) s, model: $D/output/output.cif"
        touch "$D/DONE"
    else
        echo "    build FAILED (exit $rc) after $((e - s)) s; see $D/cryoatom.log"
        echo "build exit $rc" > "$D/FAILED"
    fi
done
echo "=== batch finished $(date '+%F %T')"
