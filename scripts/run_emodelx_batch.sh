#!/bin/bash
# Run EModelX(+AF) (biomed-AI/EModelX, protocol temp_flex = map + sequence + AF templates) on every entry in
# runs/ready_entries.tsv, with its own PULCHRA and Phenix refinement options (decision A, docs/PROTOCOL.md):
#   python run.py --protocol=temp_flex --EM_map MAP --fasta FASTA --template_dir T --output_dir OUT
#                 --run_pulchra --pulchra_path modules/pulchra304/src/pulchra
#                 --run_phenix --phenix_act PHENIX --resolution RES --mul_proc_num 14
# Templates: the SAME AF3 model_0 conversions MICA used, i.e. ~/bench/runs/mica/<EMDB>/AF3_structures/<PDBID_n>/ranked_0.pdb
# (EModelX expects <template_dir>/<FASTA name>/ranked_0.pdb, chain A). --download_afdb is NOT used (it would overwrite them).
# --mul_proc_num 14 instead of the default 30: 30 forked workers would exceed the 50 GB RAM on large maps; 14 is what
# MICA chose automatically on this machine (0.6 x 24 threads). Documented deviation.
# run.py starts phenix.real_space_refine in the background (phenix.sh ends with '&'), so we wait for it.
# Resumable (DONE files), ONLY=<PDB ID> for one entry, OOM check like the MICA batch, resource gate before each entry.
#
# Usage:  bash scripts/run_emodelx_batch.sh [entries.tsv]
set -u

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
EMX=$METHODS_DIR/EModelX
RUNS=$RUNS_ROOT/emodelx
MICA_RUNS=$RUNS_ROOT/mica
LIST=${1:-$BENCH/runs/ready_entries.tsv}
TIMES=$RUNS/timings.csv
ONLY=${ONLY:-}
NPROC=${NPROC:-14}
GPU_MIN_MIB=${GPU_MIN_MIB:-12288}
RAM_MIN_GIB=${RAM_MIN_GIB:-12}

# Same threading settings as MICA (same Python 3.8 / torch 1.8.1 / MKL stack; avoids the fork deadlock and MKL error)
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_THREADING_LAYER=GNU
source "$CONDA_SH"
conda activate EModelX

mkdir -p "$RUNS"
[ -f "$TIMES" ] || echo "emdb,pdb,step,start,end,seconds,exit_code" > "$TIMES"
cd "$EMX" || exit 1

gpu_free() { nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1; }
ram_avail() { awk '/MemAvailable/ {print int($2 / 1024 / 1024)}' /proc/meminfo; }
wait_for_resources() {
    local waited=0 g r
    while true; do
        g=$(gpu_free); r=$(ram_avail)
        if [ "$g" -ge "$GPU_MIN_MIB" ] && [ "$r" -ge "$RAM_MIN_GIB" ]; then return; fi
        [ $((waited % 600)) -eq 0 ] && echo "    waiting: GPU free ${g} MiB (need ${GPU_MIN_MIB}), RAM avail ${r} GiB (need ${RAM_MIN_GIB})"
        sleep 30; waited=$((waited + 30))
    done
}

tail -n +2 "$LIST" | while IFS=$'\t' read -r N PDB CONTOUR RES RESIDUES JOBS; do
    [ -n "$N" ] || continue
    if [ -n "$ONLY" ] && [ "$PDB" != "$ONLY" ]; then continue; fi
    D=$RUNS/$N
    OUT=$D/output
    if [ -f "$D/DONE" ]; then echo "skip $PDB (done)"; continue; fi

    echo "=== $(date '+%F %T')  $PDB / EMD-$N  ($RESIDUES deposited residues)"
    rm -rf "${D:?}"
    mkdir -p "$D/templates" "$OUT"
    cp "$BENCH/dataset/entries/$N/emd_$N.map" "$D/"
    cp "$BENCH/dataset/entries/$N/${PDB,,}.fasta" "$D/"
    MISSING=0
    for H in $(grep "^>" "$D/${PDB,,}.fasta" | cut -d'|' -f1 | tr -d '>'); do
        SRC=$MICA_RUNS/$N/AF3_structures/$H/ranked_0.pdb
        if [ -s "$SRC" ]; then mkdir -p "$D/templates/$H"; cp "$SRC" "$D/templates/$H/ranked_0.pdb"; else echo "    missing template $SRC"; MISSING=1; fi
    done
    if [ "$MISSING" -eq 1 ]; then echo "templates missing" > "$D/FAILED"; continue; fi
    echo "    templates: $(ls "$D/templates" | tr '\n' ' ')"

    wait_for_resources
    OOM_BEFORE=$(awk '/^oom_kill / {print $2}' /proc/vmstat)
    s=$(date +%s)
    python run.py --protocol=temp_flex --EM_map="$D/emd_$N.map" --fasta="$D/${PDB,,}.fasta" \
        --template_dir="$D/templates" --output_dir="$OUT" \
        --run_pulchra --pulchra_path=modules/pulchra304/src/pulchra \
        --run_phenix --phenix_act="$PHENIX_ACT" --resolution="$RES" --mul_proc_num="$NPROC" \
        < /dev/null > "$D/emodelx.log" 2>&1
    rc=$?
    e=$(date +%s)
    echo "$N,$PDB,run_py,$(date -d @"$s" +%FT%T),$(date -d @"$e" +%FT%T),$((e - s)),$rc" >> "$TIMES"
    echo "    run_py: $((e - s)) s, exit $rc"
    if [ $rc -ne 0 ]; then echo "run_py exit $rc" > "$D/FAILED"; continue; fi

    # wait for the detached phenix.real_space_refine (its command line contains the output path)
    s=$(date +%s); sleep 20
    while pgrep -f "real_space_refine.*$OUT" > /dev/null; do sleep 30; done
    e=$(date +%s)
    REFINED=$(ls "$OUT"/*real_space_refined_000.* 2>/dev/null | head -1)
    echo "$N,$PDB,phenix_refine_wait,$(date -d @"$s" +%FT%T),$(date -d @"$e" +%FT%T),$((e - s)),$([ -n "$REFINED" ] && echo 0 || echo 1)" >> "$TIMES"
    echo "    phenix_refine: $((e - s)) s, refined model: ${REFINED:-MISSING}"

    OOM_AFTER=$(awk '/^oom_kill / {print $2}' /proc/vmstat)
    if [ "$OOM_AFTER" -gt "$OOM_BEFORE" ]; then
        echo "    WARNING: $((OOM_AFTER - OOM_BEFORE)) process(es) OOM-killed during this entry - results not trusted"
        echo "oom_kills $((OOM_AFTER - OOM_BEFORE))" > "$D/FAILED"
        continue
    fi
    if [ -z "$REFINED" ]; then echo "no refined model" > "$D/FAILED"; continue; fi
    touch "$D/DONE"
done
echo "=== batch finished $(date '+%F %T')"
