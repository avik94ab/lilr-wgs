#!/bin/bash
# Copy number for a whole cohort of staged slices, one SGE array task per chunk.
#
#   qsub -t 1-157 scripts/cn_array.sh config/manifest.2504.staged.tsv results/kgp2504
#
# Uses portable/lilr_cn.py rather than the Snakemake DAG because the deliverable
# is copy number: the portable script produces byte-identical calls (300/300 on
# the kgp100 cohort) without recruiting reads against the panels, which is most
# of the per-sample cost and none of what a copy number is made of.
#
# CHUNK samples per task at $NSLOTS threads. Sized for short.q: 16 samples at 4
# slots is ~10 minutes, against the 25-minute cap that keeps a job out of long.q.
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -l h_rt=00:25:00
#$ -l mem_free=4G
#$ -pe smp 4
#$ -o logs/
#$ -e logs/
#$ -N lilrcn.array
set -euo pipefail

MANIFEST=${1:?manifest}
OUTDIR=${2:?outdir}
CHUNK=${CHUNK:-16}

# `-V` carries the submitting shell's PATH, so activate the environment before
# submitting rather than naming a site path here. The module-provided samtools on
# some clusters is built without libcurl; it does not matter for staged slices,
# but a missing samtools looks like a pipeline bug three log files later.
command -v samtools >/dev/null || { echo "samtools not on PATH" >&2; exit 1; }

mkdir -p "$OUTDIR/parts"
# Header line is row 1 of the manifest, so data rows start at 2.
first=$(( (SGE_TASK_ID - 1) * CHUNK + 2 ))
last=$(( first + CHUNK - 1 ))
part="$OUTDIR/parts/chunk.$(printf '%04d' "$SGE_TASK_ID")"

# A manifest slice keeps the header, so lilr_cn.py reads it as the repository
# format rather than as a bare list of paths.
{ head -1 "$MANIFEST"; sed -n "${first},${last}p" "$MANIFEST"; } > "$part.tsv.in"
if [ "$(wc -l < "$part.tsv.in")" -le 1 ]; then
    echo "task $SGE_TASK_ID: past the end of $MANIFEST, nothing to do"
    exit 0
fi

python3 portable/lilr_cn.py -m "$part.tsv.in" -t "${NSLOTS:-4}" \
    -o "$part.tsv" --qc "$part.qc.tsv"
