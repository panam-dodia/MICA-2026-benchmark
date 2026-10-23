#!/bin/bash
# Enforce the benchmark's per-entry wall-time limit on a running batch.
# Every 5 minutes: find the entry the batch is working on (last "=== <date> <time>  <PDB> / EMD-<N>" line in its log).
# If that entry has run longer than LIMIT_H hours and has no DONE/FAILED marker, kill every process whose command line
# contains the entry's run folder (<RUNS>/<N>/) and record the reason in <RUNS>/<N>/FAILED. The batch then moves on.
# Exits when the batch log contains "=== batch finished".
#
# Usage: bash watchdog.sh <batch_log> <runs_dir> [limit_hours]
#   e.g. bash watchdog.sh ~/bench/runs/mica_batch.log ~/bench/runs/mica 12
set -u
LOG=${1:?batch log}; RUNS=${2:?runs dir}; LIMIT_H=${3:-12}
LIMIT_S=$((LIMIT_H * 3600))
echo "$(date '+%F %T') watchdog started: $LOG, $RUNS, limit ${LIMIT_H} h"

while true; do
    if tail -1 "$LOG" 2>/dev/null | grep -q "^=== batch finished"; then
        echo "$(date '+%F %T') batch finished; watchdog exiting"; exit 0
    fi
    LAST=$(grep -E "^=== [0-9]{4}-[0-9]{2}-[0-9]{2} " "$LOG" | tail -1)
    if [ -n "$LAST" ]; then
        START=$(echo "$LAST" | awk '{print $2" "$3}')
        PDB=$(echo "$LAST" | awk '{print $4}')
        N=$(echo "$LAST" | sed -E 's/.*EMD-([0-9]+).*/\1/')
        D=$RUNS/$N
        ELAPSED=$(( $(date +%s) - $(date -d "$START" +%s) ))
        if [ "$ELAPSED" -gt "$LIMIT_S" ] && [ ! -f "$D/DONE" ] && [ ! -f "$D/FAILED" ]; then
            echo "$(date '+%F %T') $PDB (EMD-$N) exceeded ${LIMIT_H} h (${ELAPSED} s); stopping it"
            pkill -f "$D/"
            sleep 15
            echo "timeout: exceeded ${LIMIT_H} h per-entry limit (${ELAPSED} s since $START); stopped by watchdog $(date '+%F %T')" > "$D/FAILED"
        fi
    fi
    sleep 300
done
