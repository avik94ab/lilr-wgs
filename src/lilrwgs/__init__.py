"""LILR genotyping from short-read whole-genome sequencing.

Stages are independent modules with thin CLI wrappers, so each can be run and
tested without Snakemake and debugged without a cluster:

    loci         GRCh38 coordinates for the leukocyte receptor complex
    extract      CRAM -> LRC read pairs -> FASTQ
    coverage     per-sample depth model: lambda_1, dispersion, GC, ALT verdict
    depth_model  callability thresholds derived from that model
    cn           copy number per locus
    assign       reads -> genes, arbitrating the paralogues
    genotype     variants, phasing, consensus
    sequences    CDS / cDNA / protein from a genomic haplotype
"""

__version__ = "0.1.0"
