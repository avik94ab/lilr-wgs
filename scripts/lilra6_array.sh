#!/bin/bash
# LILRA6 copy number by extract-and-realign, one SGE array task per chunk.
#
#   qsub -t 1-25 scripts/lilra6_array.sh inputs.txt results/lilra6
#
# Chunked rather than one task per sample because the per-task cost is dominated
# by loading the bwa index, not by aligning: the index is 5.3 GB and the read set
# is a few hundred thousand reads. One sample per task would pay that load 100
# times; CHUNK samples per task pays it CHUNK times fewer -- bwa is invoked once
# per sample either way, but a warm page cache makes every load after the first
# nearly free.
#
# Sized for short.q at CHUNK=4: measured ~2 min per sample after the first, so a
# task is ~10 minutes against the 25-minute cap that keeps a job out of long.q.
#
# mem_free on SGE is *per slot*, so 2G at 8 slots requests 16 GB -- comfortably
# above the ~5.3 GB the index occupies in one bwa process. Asking for 10G here
# would request 80 GB and queue behind nothing but its own arithmetic.
#
# --jobs is deliberately left at 1 inside a task: two bwa processes in one task
# means two copies of a 5.3 GB index, which is a worse use of the same slots
# than two tasks.
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -l h_rt=00:25:00
#$ -l mem_free=2G
#$ -pe smp 8
#$ -o logs/
#$ -e logs/
#$ -N lilra6.array
set -euo pipefail

INPUTS=${1:?inputs list}
OUTDIR=${2:?outdir}
REFERENCE=${REFERENCE:-resources/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa}
# What to realign to. Defaults to REFERENCE, i.e. GRCh38. Set TARGET to
# resources/reference/chm13v2.0.fa to measure in T2T coordinates, where LILRA3
# is on the primary assembly rather than on four alt contigs.
TARGET=${TARGET:-$REFERENCE}
# Which genes to report. Empty means everything the target can do.
GENES=${GENES:-}
CHUNK=${CHUNK:-4}

command -v bwa >/dev/null || { echo "bwa not on PATH" >&2; exit 1; }
command -v samtools >/dev/null || { echo "samtools not on PATH" >&2; exit 1; }
[ -s "$TARGET.bwt" ] || { echo "$TARGET.bwt missing; build or fetch the bwa index" >&2; exit 1; }
# The failure with no error message. bwa does not report a missing .alt file; it
# aligns without ALT-awareness and every LILRA6 call comes back not_measured, so
# a whole array would burn its slots producing refusals. Only GRCh38 has ALT
# contigs -- on CHM13 the file's absence is correct, not a misconfiguration.
#
# Which assembly this is comes from chr19's length in the .fai, not from the
# filename: renaming chm13v2.0.fa must not be able to change what it is measured
# as, and lilrwgs.realign.assembly_of_reference decides it the same way.
[ -s "$TARGET.fai" ] || { echo "$TARGET.fai missing; samtools faidx it" >&2; exit 1; }
CHR19_LEN=$(awk -F'\t' '$1=="chr19"{print $2}' "$TARGET.fai")
if [ "$CHR19_LEN" = "58617616" ]; then
    [ -s "$TARGET.alt" ] || { echo "$TARGET.alt missing; run scripts/fetch_bwa_index.sh" >&2; exit 1; }
fi

# Retried, because `mkdir -p` is not reliably idempotent across nodes on a
# parallel filesystem and neither is the `-d` test that would check it. 25 array
# tasks starting at once all raced to create these: on the first run three lost
# with "File exists", and under `set -e` that killed the task before it did any
# work, so the cohort came back 88/100 with three chunks simply absent. Merely
# tolerating EEXIST was not enough either -- the next run still lost one task,
# this time to a `-d` that returned false for a directory another node had
# already created, because BeeGFS had not propagated the metadata yet.
#
# So: try, sleep, look again. A directory that is genuinely uncreatable (no
# permission, full quota) still fails loudly after the last attempt.
for attempt in 1 2 3 4 5; do
    mkdir -p "$OUTDIR/parts" logs 2>/dev/null || true
    missing=""
    for d in "$OUTDIR/parts" logs; do
        [ -d "$d" ] || missing="$missing $d"
    done
    [ -z "$missing" ] && break
    if [ "$attempt" -eq 5 ]; then
        echo "cannot create:$missing" >&2
        exit 1
    fi
    sleep $(( attempt * 2 ))
done

first=$(( (SGE_TASK_ID - 1) * CHUNK + 1 ))
last=$(( first + CHUNK - 1 ))
part="$OUTDIR/parts/chunk.$(printf '%04d' "$SGE_TASK_ID")"

sed -n "${first},${last}p" "$INPUTS" > "$part.inputs.txt"
if [ ! -s "$part.inputs.txt" ]; then
    echo "task $SGE_TASK_ID: past the end of $INPUTS, nothing to do"
    exit 0
fi

# $PYTHONPATH unset is the normal state on a compute node and fatal under
# `nounset`, which this script runs with.
export PYTHONPATH="src${PYTHONPATH:+:$PYTHONPATH}"

python3 scripts/lilra6_cn.py \
    --inputs "$part.inputs.txt" \
    --reference "$REFERENCE" \
    --target "$TARGET" \
    --threads "${NSLOTS:-8}" \
    --jobs 1 \
    ${GENES:+--genes "$GENES"} \
    --outdir "$OUTDIR/qc/$(printf '%04d' "$SGE_TASK_ID")" \
    -o "$part.tsv"
