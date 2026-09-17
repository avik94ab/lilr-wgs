#!/bin/bash
# Fetch the GRCh38 the 1000 Genomes CRAMs were compressed against.
#
# It must be *this* reference. A different GRCh38 build decodes most bases
# correctly and corrupts the rest, which is considerably worse than failing --
# the CRAM's @SQ M5 checksums are what make the difference detectable, and
# samtools will not always check them.
set -euo pipefail
DEST="${1:-resources/reference}"
URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
mkdir -p "$DEST"
FA="$DEST/GRCh38_full_analysis_set_plus_decoy_hla.fa"

if [ -s "$FA" ]; then
    echo "$FA present ($(du -h "$FA" | cut -f1))"
else
    echo "fetching ~3.2 GB from EBI ..."
    curl -fL --retry 3 -o "$FA" "$URL"
fi
[ -s "$FA.fai" ] || samtools faidx "$FA"
echo "reference ready: $FA"
