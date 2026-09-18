#!/bin/bash
# Assert the environment can actually do the job, with reasons.
#
# The one that matters is libcurl. Without it `samtools view https://...` fails
# with a bare "fail to open file", which reads exactly like a typo in a path and
# sends people looking in the wrong place for an afternoon.
set -uo pipefail
fail=0

for tool in samtools bcftools bwa bowtie2 bowtie2-build gatk whatshap miniprot python; do
    if command -v "$tool" > /dev/null; then
        printf "  %-14s %s\n" "$tool" "$(command -v "$tool")"
    else
        printf "  %-14s MISSING\n" "$tool"; fail=1
    fi
done

if command -v samtools > /dev/null; then
    # samtools prints a Features line for itself and another for htslib;
    # libcurl belongs to htslib, so grep the whole output, not the first line.
    if samtools --version | grep -q "libcurl=yes"; then
        echo "  samtools has libcurl: remote CRAM reads will work"
    else
        echo "  ERROR: this samtools is built without libcurl."
        echo "         It cannot read CRAMs over HTTPS, which is the front end of"
        echo "         this pipeline. Use environment.yml rather than a system or"
        echo "         module build, or stage CRAMs locally."
        fail=1
    fi
fi

# The `.alt` file, which bwa finds by name and never mentions. Without it the
# realignment is not ALT-aware, every MAPQ-20 window in the LRC reads near zero
# for every sample alike, and a whole cohort comes out as apparent LILRA6
# deletion homozygotes. There is no error message from bwa for this -- the run
# simply succeeds and is wrong -- so it is checked here instead.
REF="${LILRWGS_REFERENCE:-resources/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa}"
if [ -s "$REF.alt" ]; then
    echo "  bwa .alt present: realignment will be ALT-aware"
elif [ -s "$REF" ]; then
    echo "  WARNING: $REF.alt is missing."
    echo "           bwa gives no error for this; it just aligns without"
    echo "           ALT-awareness, and every LILRA6 call becomes not_measured."
    echo "           Run scripts/fetch_bwa_index.sh"
fi

if [ "$fail" -ne 0 ]; then
    echo; echo "Environment is incomplete. micromamba create -y -f environment.yml"
    exit 1
fi
echo; echo "Environment OK."
