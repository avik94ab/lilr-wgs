#!/bin/bash
# Report coherence against a baseline as an array job's chunks land.
#
#   bash scripts/watch_coherence.sh <parts-dir> <baseline-tsv> <gene> [job-id]
#
# Emits one line per reporting boundary, so a systematic problem surfaces while
# there is still time to kill the run rather than after it. Reports densely at
# the start -- every chunk for the first five -- and then every fifth chunk,
# because the failure modes this is watching for are systematic: if the first
# hundred samples agree, the next two thousand are not going to start
# disagreeing for a new reason, and a report every two minutes for eight hours
# is noise that gets itself ignored.
set -uo pipefail

PARTS=${1:?parts dir}
BASELINE=${2:?baseline tsv}
GENE=${3:-LILRA6}
JOB=${4:-}

DENSE_UNTIL=5     # report on every chunk up to this many
THEN_EVERY=5      # and every Nth chunk after

last=0
while true; do
    n=$(ls "$PARTS"/chunk.[0-9]*.tsv 2>/dev/null | wc -l)

    if [ "$n" -gt "$last" ]; then
        report=0
        if [ "$n" -le "$DENSE_UNTIL" ]; then
            report=1
        elif [ $(( n % THEN_EVERY )) -eq 0 ]; then
            report=1
        fi
        if [ "$report" = "1" ]; then
            merged=$(mktemp)
            head -1 "$(ls "$PARTS"/chunk.[0-9]*.tsv | head -1)" > "$merged"
            tail -qn +2 "$PARTS"/chunk.[0-9]*.tsv >> "$merged" 2>/dev/null
            python3 validation/compare_realign.py --realigned "$merged" \
                --as-is "$BASELINE" --gene "$GENE" 2>/dev/null \
              | awk -v n="$n" -v g="$GENE" '
                  /integer agreement/ {agree=$3; pct=$4}
                  /^  median/ {med=$2}
                  /both measured/ {meas=$3}
                  END {printf "[%s chunks] %s: %s %s  measured %s  median delta %s\n",
                               n, g, agree, pct, meas, med}'
            rm -f "$merged"
        fi
        last=$n
    fi

    # Stop when the job is gone, after one final report.
    if [ -n "$JOB" ] && ! qstat -j "$JOB" > /dev/null 2>&1; then
        echo "[done] $JOB has left the queue with $n chunks"
        exit 0
    fi
    sleep 60
done
