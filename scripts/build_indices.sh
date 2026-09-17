#!/bin/bash
# build_indices.sh PANEL_DIR PANEL_IDX LOCUS_REFS LOCUS_IDX THREADS
#
# bowtie2 indices for the 11 HPRC pangenome panels (read recruitment) and the 11
# single per-locus references (variant calling). Built rather than shipped: they
# are ~30 MB of binary blobs derived from FASTAs that are already in the repo.
#
# Two reference sets, deliberately. Recruitment uses a 465-sequence pangenome
# panel so that a divergent allele still attracts its own reads; calling uses one
# reference per locus so that variants have a single coordinate system. Using the
# panel for calling would give no usable VCF, and using the single reference for
# recruitment loses reads from every haplotype unlike it.
set -euo pipefail
PANEL_DIR="${1:?panel dir}"; PANEL_IDX="${2:?panel index dir}"
LOCUS_REFS="${3:?locus ref dir}"; LOCUS_IDX="${4:?locus index dir}"
THREADS="${5:-4}"

mkdir -p "$PANEL_IDX" "$LOCUS_IDX"

for fa in "$PANEL_DIR"/*.fasta; do
    gene=$(basename "$fa" .fasta)
    if [ -f "$PANEL_IDX/${gene}.1.bt2" ]; then
        echo "panel $gene: index present, skipping"; continue
    fi
    echo "panel $gene: building"
    bowtie2-build --threads "$THREADS" "$fa" "$PANEL_IDX/$gene" > /dev/null
done

for fa in "$LOCUS_REFS"/*_named.fa; do
    gene=$(basename "$fa" _named.fa)
    if [ -f "$LOCUS_IDX/${gene}.1.bt2" ]; then
        echo "locus $gene: index present, skipping"; continue
    fi
    echo "locus $gene: building"
    bowtie2-build --threads "$THREADS" "$fa" "$LOCUS_IDX/$gene" > /dev/null
    samtools faidx "$fa"
done

echo "indices ready"
