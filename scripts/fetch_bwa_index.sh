#!/bin/bash
# Fetch the bwa index for the GRCh38 analysis set, including the `.alt` file.
#
# Downloaded rather than built, for a reason that is not about the ~90 minutes
# `bwa index` would take. EBI publishes the index NYGC aligned the 1000 Genomes
# 30x CRAMs against, so realigning against *this* index reproduces their
# placements rather than merely resembling them -- which is what makes the
# realigned copy numbers comparable to the ones called from the CRAMs as-is.
#
# The `.alt` file is the load-bearing one. `bwa mem` looks for `{idxbase}.alt`
# and, finding it, computes MAPQ over primary hits alone. Without it every read
# in the LRC has nine more equally good placements on the alt haplotypes and
# comes back at MAPQ 0 -- and every MAPQ-20 window this pipeline measures reads
# near zero, for everyone, which looks exactly like a cohort of homozygous
# deletions. scripts/check_env.sh asserts the file is present.
set -euo pipefail
DEST="${1:-resources/reference}"
BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome"
FA="$DEST/GRCh38_full_analysis_set_plus_decoy_hla.fa"

if [ ! -s "$FA" ]; then
    echo "error: $FA not found; run scripts/fetch_reference.sh first" >&2
    exit 1
fi

mkdir -p "$DEST"
# .alt first: it is 476 KB against 5.3 GB for the rest, so a run interrupted
# part-way through still leaves the file whose absence is silently wrong rather
# than loudly missing.
for ext in alt amb ann pac bwt sa; do
    target="$FA.$ext"
    if [ -s "$target" ]; then
        echo "  $(basename "$target") present ($(du -h "$target" | cut -f1))"
        continue
    fi
    echo "  fetching $(basename "$target") ..."
    curl -fL --retry 3 --retry-delay 5 -o "$target.part" \
        "$BASE/GRCh38_full_analysis_set_plus_decoy_hla.fa.$ext"
    # Rename only on success: a truncated .bwt is accepted by bwa's loader far
    # enough to produce garbage alignments, and a partial file left in place
    # would be picked up by the next run as complete.
    mv "$target.part" "$target"
done

echo
echo "bwa index ready: $FA"
echo "  .alt present -> bwa mem will align ALT-aware"
