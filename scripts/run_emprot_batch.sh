#!/bin/bash
# Run EMProt v1.2 (map + sequence + AF3 single-chain templates) on every entry in runs/ready_entries.tsv:
#   emprot build --map MAP --seq seq.fa --chain t1.cif t2.cif ... --device 0 --output OUT
# Inputs from scripts/13_prepare_emprot_inputs.py: one FASTA record and one AF3 model_0 template per chain copy (the
# same AF3 files MICA and EModelX(+AF) use). Templates are copied to the run folder on ext4, as EMProt uses symlinks.
# Resumable: entries with a DONE file are skipped. Optional: ONLY=<PDB ID> runs one entry.
# Before each entry it waits for >= GPU_MIN_MIB free GPU memory and >= RAM_MIN_GIB available RAM (safety check).
# Final model: <entry>/output/output.cif (unrefined, EMProt's default output).
#
# Usage:  bash scripts/run_emprot_batch.sh [entries.tsv]
set -u

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
RUNS=$RUNS_ROOT/emprot
LIST=${1:-$BENCH/runs/ready_entries.tsv}
TIMES=$RUNS/timings.csv
ONLY=${ONLY:-}
GPU_MIN_MIB=${GPU_MIN_MIB:-12288}
RAM_MIN_GIB=${RAM_MIN_GIB:-12}

source "$CONDA_SH"
conda activate emprot
# EMProt README: avoids GLIBCXX errors from the system libstdc++
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}

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
    cp "$BENCH/runs/emprot_inputs/$N/seq.fa" "$D/"
    TEMPLATES=()
    k=0
    while read -r T; do
        [ -n "$T" ] || continue
        k=$((k + 1))
        cp "$BENCH/$T" "$D/template_$k.cif"
        TEMPLATES+=("$D/template_$k.cif")
    done < "$BENCH/runs/emprot_inputs/$N/templates.txt"
    echo "    chain copies / templates: $k"

    wait_for_resources
    s=$(date +%s)
    emprot build --map "$D/emd_$N.map" --seq "$D/seq.fa" --chain "${TEMPLATES[@]}" --device 0 --output "$D/output" < /dev/null > "$D/emprot.log" 2>&1
    rc=$?
    e=$(date +%s)
    echo "$N,$PDB,build,$(date -d @"$s" +%FT%T),$(date -d @"$e" +%FT%T),$((e - s)),$rc" >> "$TIMES"
    if [ $rc -eq 0 ] && [ -s "$D/output/output.cif" ]; then
        echo "    build: $((e - s)) s, model: $D/output/output.cif"
        touch "$D/DONE"
    else
        echo "    build FAILED (exit $rc) after $((e - s)) s; see $D/emprot.log"
        echo "build exit $rc" > "$D/FAILED"
    fi
done
echo "=== batch finished $(date '+%F %T')"
