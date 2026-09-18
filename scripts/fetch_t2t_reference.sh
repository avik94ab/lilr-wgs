#!/bin/bash
# Fetch T2T-CHM13v2.0 and build its bwa index.
#
# Why a second reference at all: GRCh38's chr19 carries the common ~6.7 kb
# LILRA3 deletion, so LILRA3 is annotated only on alt contigs and has to be read
# as MAPQ-0 depth over four near-identical haplotypes, counting supplementary
# records. That route cannot survive a regional extraction -- the reads it counts
# have primaries scattered genome-wide -- so the realigning front end cannot call
# LILRA3 at all (PLAN.md S12).
#
# CHM13v2.0 carries the insertion allele. Measured here rather than assumed: our
# 7,126 bp LILRA3 calling reference aligns to chr19:57,377,084-57,384,209 at
# 7,126/7,126 identity, reverse-complemented, mapq 60, with every other hit in
# the region a partial paralogue at 59-87%. So on CHM13 LILRA3 is ordinary
# single-copy primary sequence, measurable at MAPQ 20 like any other gene.
#
# Note the annotation does NOT list LILRA3. chm13v2.0_RefSeq_Liftoff_v5.1.gff3
# is lifted from the GRCh38 primary assembly, which has no LILRA3 to lift, so its
# absence there says nothing about the sequence. Do not take that file as
# evidence either way.
#
# Unlike GRCh38, no prebuilt bwa index is published for CHM13 (checked: the
# analysis_set and indexes/ paths 404), so it is built here -- roughly 60-90
# minutes and ~4.5 GB. There is no `.alt` file and there should not be: CHM13 has
# no ALT contigs, which is most of the point.
set -euo pipefail
DEST="${1:-resources/reference}"
BASE="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies"
FA="$DEST/chm13v2.0.fa"

mkdir -p "$DEST"

if [ -s "$FA" ]; then
    echo "$FA present ($(du -h "$FA" | cut -f1))"
else
    echo "fetching ~936 MB (gzipped) from the human-pangenomics bucket ..."
    curl -fL --retry 3 -o "$FA.gz.part" "$BASE/analysis_set/chm13v2.0.fa.gz"
    mv "$FA.gz.part" "$FA.gz"
    echo "decompressing ..."
    gunzip -c "$FA.gz" > "$FA.part"
    mv "$FA.part" "$FA"
fi

[ -s "$FA.fai" ] || samtools faidx "$FA"

# The annotation, for deriving coordinates. Small, and the provenance for every
# interval in lilrwgs.loci_chm13 except LILRA3's.
GFF="$DEST/chm13v2.0_RefSeq_Liftoff_v5.1.gff3.gz"
if [ ! -s "$GFF" ]; then
    echo "fetching the annotation ..."
    curl -fL --retry 3 -o "$GFF.part" \
        "$BASE/annotation/chm13v2.0_RefSeq_Liftoff_v5.1.gff3.gz"
    mv "$GFF.part" "$GFF"
fi

if [ -s "$FA.bwt" ]; then
    echo "bwa index present"
else
    echo "building the bwa index (60-90 min, ~4.5 GB) ..."
    bwa index "$FA"
fi

echo
echo "T2T reference ready: $FA"
echo "  no .alt file, and none is wanted: CHM13 has no ALT contigs"
