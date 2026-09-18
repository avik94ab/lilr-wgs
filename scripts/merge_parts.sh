#!/bin/bash
# Concatenate an array job's per-chunk TSVs into one, keeping a single header.
#
#   bash scripts/merge_parts.sh results/lilra6/parts results/lilra6/lilra6_cn.tsv
#
# Counts what it merged and says so, because the interesting failure is a chunk
# that produced no file at all -- a task killed by the scheduler leaves nothing
# behind, and a merged table that is quietly short by four samples looks exactly
# like a merged table that is complete.
set -euo pipefail
PARTS=${1:?parts dir}
OUT=${2:?output tsv}
PATTERN=${3:-chunk.*.tsv}

shopt -s nullglob
files=("$PARTS"/$PATTERN)
if [ ${#files[@]} -eq 0 ]; then
    echo "no files matching $PARTS/$PATTERN" >&2
    exit 1
fi

mkdir -p "$(dirname "$OUT")"
head -1 "${files[0]}" > "$OUT"
for f in "${files[@]}"; do
    tail -n +2 "$f"
done >> "$OUT"

echo "merged ${#files[@]} parts -> $OUT ($(( $(wc -l < "$OUT") - 1 )) rows)"
