# Validation

The HPRC panels this pipeline aligns against carry their donor in the FASTA
header, so the same files are a truth set: **232 donors, 101 of whom also have a
1000 Genomes 30× CRAM**.

```bash
python validation/build_truth.py --panels resources/gdna \
    --kgp-samples validation/kgp_2504_samples.txt --outdir validation/truth
python scripts/make_manifest.py --collection 2504 \
    --samples validation/truth/overlap_1kgp.txt -o config/manifest.hprc.tsv

# 1. Stage the slices. On a cluster whose compute nodes have no outbound route —
#    Wynton's do not — this has to happen somewhere that does. Skip it and the
#    jobs fail with "Destination address required" rather than with anything
#    that mentions the network.
python3 scripts/stage_slices.py --manifest config/manifest.hprc.tsv \
    --reference resources/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    --outdir $SCRATCH/lilr_slices \
    --out-manifest config/manifest.hprc.staged.tsv --jobs 8

# 2. Copy number for all 101. `cohort_scale.json` is the target rather than
#    `summary.csv` because the genotype rules are not needed to score CN.
snakemake -s workflow/Snakefile --profile profiles/sge_staged \
    --config manifest=config/manifest.hprc.staged.tsv outdir=results/asis \
    --jobs 40 \
    $PWD/results/asis/cn/cohort_scale.json $PWD/results/asis/qc/coverage_cohort.tsv

python validation/compare_cn.py results/asis/cn/cn_calls.tsv \
    -o validation/reports/overlap101_asis.txt
```

Two things about that Snakemake line that cost an hour between them, so they are
written out rather than left to be rediscovered:

* **The targets must come after a flag that is not `--config`.** `--config` takes
  a variable number of arguments, so targets placed directly after it are parsed
  as config entries and the run dies on `Invalid config definition`. `--jobs 40`
  sits between them here for exactly that reason.
* **The targets must be absolute paths.** `outdir` resolves to an absolute path
  inside the workflow, and a relative target does not match it — the failure is
  `MissingRuleException`, which points at the rules rather than at the path.

## Read the caveats, they are not decoration

**Circularity.** A donor in the overlap is aligned against a panel containing its
own haplotypes, which inflates recruitment and accuracy relative to an unseen
sample. `build_truth.py --leave-one-donor-out <dir>` writes donor-excluded
panels; rerunning against those and scoring with `--leave-one-donor-out` gives
the number that generalises. Report both, labelled.

```bash
python validation/build_truth.py --panels resources/gdna \
    --kgp-samples validation/kgp_2504_samples.txt --outdir validation/truth \
    --leave-one-donor-out $SCRATCH/lilr_lodo/panels

# `{sample}` in either panel path makes recruitment per-donor, and turns
# build_indices into a per-donor rule. The slices are the same ones staged above.
snakemake -s workflow/Snakefile --profile profiles/sge_staged \
    --config manifest=config/manifest.hprc.staged.tsv outdir=results/lodo \
             gdna_panels=$SCRATCH/lilr_lodo/panels/{sample} \
             panel_index=$SCRATCH/lilr_lodo/index/{sample} \
    --jobs 40 \
    $PWD/results/lodo/cn/cohort_scale.json $PWD/results/lodo/qc/coverage_cohort.tsv

python validation/compare_cn.py results/lodo/cn/cn_calls.tsv \
    --leave-one-donor-out -o validation/reports/overlap101_lodo.txt
```

A per-donor run builds recruitment panels only — the per-locus calling
references do not depend on the donor, and 101 jobs racing to build one shared
index is a way to get a half-written one. So `scripts/build_indices.sh` must have
been run once beforehand; the workflow refuses to start if it has not been.

**Inferred absences.** LILRA3 CN 0 in the truth set comes from a donor being
absent from that panel while present in at least 8 of the other 10, not from
counting copies. Score with `--counted-only` to exclude those rows.

**Confident vs flagged.** `compare_cn.py` scores them separately on purpose. A
method that is 99% accurate on its confident calls and 60% on its flagged ones is
more useful than one that is uniformly 95%, and the difference is invisible in a
single accuracy number. But note that a flag firing unevenly across copy-number
classes biases any allele frequency computed from the confident subset only —
compute frequencies on the unfiltered calls.

## Reports

- [`reports/pilot5_cn.txt`](reports/pilot5_cn.txt) — five donors spanning LILRA6
  CN 1–4 and LILRA3 CN 0/2, 15/15 correct. The run that found the LILRA3
  supplementary-alignment bug: before the truth set existed, a uniform
  "LILRA3 CN 0" across every sample looked entirely plausible.
