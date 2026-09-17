# lilr-wgs

Copy number and phased allele sequences for the 11 **LILR** genes, called from
**short-read whole-genome sequencing** — specifically the 1000 Genomes Project
30× GRCh38 CRAMs.

It is a sibling of [`lilr-genotyper`](https://github.com/avik94ab/lilr-genotyper),
which does the same job from targeted-capture FASTQs at several hundred ×. Almost
every depth-related constant in that pipeline is calibrated to that coverage, and
at 30× those constants stop being filters and become masks: a `DP>=20` variant
filter discards the majority of true heterozygous sites, and a `depth < 15 → N`
consensus mask blanks roughly half of every gene.

So the organising idea here is that **read depth is the subject, not a
parameter**. Every threshold is derived from the sample's own measured coverage,
is aware of the locus's copy number, and is two-sided — because in a cluster of
~90%-identical paralogues, a position at three times the expected depth is a
pile-up of reads that the assignment step failed to separate, and a variant
called there is a paralogous sequence variant wearing a heterozygote's clothes.

See [`PLAN.md`](PLAN.md) for the design and the phase-by-phase build, and
[`docs/`](docs/) for the method notes.

## Status

Under construction. `PLAN.md` tracks what is done and what is not; nothing here
should be treated as validated until Phase 8 reports numbers.

## Why srWGS changes the problem

| | capture (`lilr-genotyper`) | srWGS 30× (here) |
|---|---|---|
| input | whole-readset FASTQ | GRCh38 CRAM, aligned, ALT-aware |
| front end | 11 full-readset bowtie2 passes | one CRAM slice over the LRC |
| duplicates | must **not** be marked — probe geometry, not PCR | already marked; honoured |
| depth baseline | anchor LILR loci inside the same panel | control loci outside the cluster — a genuinely external baseline |
| copy number | cohort GMM, needs ≥30 samples and a DAG barrier | absolute and per-sample; the cohort pass is a refinement |
| depth thresholds | fixed integers | quantiles of a fitted per-copy depth distribution |

The baseline is the deep change. With capture, the only available reference for
"what does one copy look like?" is other LILR genes in the same panel, so copy
number is inherently relative and a cohort has to be called as one batch. With
WGS the answer is measurable from each sample's own genome, so a single sample
can be called on its own.

## The region is difficult, and the difficulties are specific

Three things about the leukocyte receptor complex on chr19q13.42 shape every
design decision, and the code says so where it matters:

1. **LILRA3 is not in the GRCh38 primary assembly.** The reference chromosome
   carries the common ~6.7 kb deletion, so LILRA3 is annotated only on LRC alt
   contigs. Depth over LILRA3 means depth over those contigs at MAPQ 0, plus an
   independent assay at the deletion junction, chr19:54,296,977.
2. **LILRA6 and LILRB3 are ~97% identical where the reads are.** Only the 3'
   ends are paralogue-unique — 2,900 bp and 1,881 bp. Short reads cannot be
   attributed to one gene in the 5' block, so they are kept for both rather than
   competitively discarded, which is what PING does for KIR2DL5A/B.
3. **Nine alt haplotypes are all placed at one primary interval** that contains
   every window used here. ALT-aware alignment makes those windows measurable;
   non-ALT-aware alignment makes every MAPQ-20 count read near zero, which looks
   exactly like a homozygous deletion. The pipeline measures which case it is
   per sample instead of trusting the specification.

## Requirements

```bash
micromamba create -y -f environment.yml && micromamba activate lilr-wgs
bash scripts/check_env.sh          # asserts samtools has libcurl, etc.
```

`samtools` **must** be built with libcurl — reading CRAMs over HTTPS is the
front end. The check script fails loudly if it is not, because without it the
failure surfaces as a bare "fail to open file" that reads like a bad path.

## Provenance

- Cross-mapping arbitration, the HPRC gene panels and the sequence builder come
  from `lilr-genotyper`.
- The GRCh38 coordinates, the LILRA3 junction assay and the ALT-awareness
  diagnostic come from `lilrCN_aou`.
- The extractor pattern, ratio-to-a-reference-locus copy number with a
  human-overridable threshold file, and the precedent for not competitively
  filtering an inseparable paralogue pair come from
  [PING](https://github.com/Hollenbach-lab/PING).
