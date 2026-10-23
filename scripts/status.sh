#!/bin/bash
# One-screen status of the benchmark. Usage: bash scripts/status.sh
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
R=$RUNS_ROOT
RAW=$BENCH/eval/raw
echo "=== $(date '+%F %T')   uptime since $(uptime -s)"

echo; echo "### RUNNING / QUEUED"
P=$(pgrep -af "run_.*_batch.sh|while pgrep|watchdog.sh|score_raw.sh" | grep -v pgrep | awk '{$1=""; print}' | sort | uniq -c | cut -c1-120)
[ -n "$P" ] && echo "$P" || echo "  nothing running"

echo; echo "### METHODS   (done / failed / total 30)"
for m in mica modelangelo cryoatom emprot emodelx; do
    d=$(ls $R/$m/*/DONE 2>/dev/null | wc -l)
    f=$(ls $R/$m/*/FAILED 2>/dev/null | wc -l)
    cur=""
    if [ -f $R/${m}_batch.log ]; then
        last=$(grep -E '^===' $R/${m}_batch.log | tail -1)
        case "$last" in *"batch finished"*) cur="finished ${last#=== batch finished }";; *) cur="now: $(echo "$last" | cut -c5-45)";; esac
    fi
    printf "  %-12s done=%-3s failed=%-3s %s\n" $m $d $f "$cur"
done

echo; echo "### FAILURES (reason)   [stale] = left from an earlier attempt; the current batch has not re-run it yet"
for m in mica modelangelo cryoatom emprot emodelx; do
    START=$(grep -m1 -oE "^=== [0-9-]+ [0-9:]+" $R/${m}_batch.log 2>/dev/null | cut -c5-)
    # no current batch log yet -> the method has not started its current attempt, so every marker is stale
    START_S=$( [ -n "$START" ] && date -d "$START" +%s || date +%s )
    for f in $R/$m/*/FAILED; do
        [ -f "$f" ] || continue
        tag=""; [ "$(stat -c %Y "$f")" -lt "$START_S" ] && tag="[stale] "
        echo "  $m $(basename $(dirname $f)): ${tag}$(head -c 90 $f)"
    done
done

echo; echo "### SCORING   (models scored)"
for m in mica modelangelo cryoatom emprot emodelx_af; do
    n=$(ls $RAW/$m/*.cc_recall.txt 2>/dev/null | wc -l); u=$(grep -l "^TM-score=" $RAW/$m/*.usalign.txt 2>/dev/null | wc -l)
    printf "  %-12s chain_comparison=%-3s US-align=%-3s\n" $m $n $u
done

echo; echo "### MACHINE"
echo "  GPU: $(nvidia-smi --query-gpu=memory.used,memory.total,utilization.gpu --format=csv,noheader)"
echo "  RAM: $(free -g | awk '/Mem:/ {print $3" GB used of "$2" GB, available "$7" GB"}')"
echo "  OOM kills since boot: $(awk '/^oom_kill / {print $2}' /proc/vmstat)"
