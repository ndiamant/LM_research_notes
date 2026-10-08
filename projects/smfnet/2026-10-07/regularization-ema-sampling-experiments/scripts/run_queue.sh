#!/usr/bin/env bash
# Run queued overnight runs one at a time. Each line of QUEUE is
# "NAME override ...". The next line is read only when a run finishes, so the
# queue can be edited while runs are going. Waits first for any run listed as
# started but not ended in runs.txt (e.g. one launched by hand), and stops
# taking new runs once fewer than MIN_LEFT_MIN minutes remain in the job.
set -uo pipefail

ROOT="${SCRATCH}/smf_models/overnight_20261006"
QUEUE="${QUEUE:-$ROOT/queue.txt}"
MIN_LEFT_MIN="${MIN_LEFT_MIN:-50}"
HERE="$(cd "$(dirname "$0")" && pwd)"

running() {
  local started ended
  started=$(grep -c " start " "$ROOT/logs/runs.txt" 2>/dev/null || true)
  ended=$(grep -c " end " "$ROOT/logs/runs.txt" 2>/dev/null || true)
  [ "${started:-0}" -gt "${ended:-0}" ]
}

minutes_left() {
  local end
  end=$(scontrol show job "$SLURM_JOB_ID" | grep -oP 'EndTime=\K\S+')
  echo $(( ($(date -d "$end" +%s) - $(date +%s)) / 60 ))
}

while running; do sleep 120; done
while true; do
  line=$(grep -v '^\s*\(#\|$\)' "$QUEUE" | head -1)
  [ -z "$line" ] && { echo "[$(date '+%F %T')] queue empty" >> "$ROOT/logs/runs.txt"; break; }
  left=$(minutes_left)
  if [ "$left" -lt "$MIN_LEFT_MIN" ]; then
    echo "[$(date '+%F %T')] stopping: ${left} min left in job" >> "$ROOT/logs/runs.txt"
    break
  fi
  # Remove the line before running it, so a crash does not repeat it.
  grep -vxF "$line" "$QUEUE" > "$QUEUE.tmp" && mv "$QUEUE.tmp" "$QUEUE"
  # shellcheck disable=SC2086
  "$HERE/run_one.sh" $line
done
