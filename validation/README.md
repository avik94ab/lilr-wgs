# Validation

The HPRC panels this pipeline aligns against carry their donor in the FASTA
header, so the same files are a truth set: **232 donors, 101 of whom also have a
1000 Genomes 30× CRAM**.

```bash
python validation/build_truth.py --panels resources/gdna \
    --kgp-samples validation/kgp_2504_samples.txt --outdir validation/truth
python scripts/make_manifest.py --collection 2504 \
    --samples validation/truth/overlap_1kgp.txt -o config/manifest.hprc.tsv

snakemake -s workflow/Snakefile --configfile config/config.yaml \
    --config manifest=config/manifest.hprc.tsv --profile profiles/sge

python validation/compare_cn.py results/cn/cn_calls.tsv -o validation/reports/overlap101.txt
```

## Read the caveats, they are not decoration

**Circularity.** A donor in the overlap is aligned against a panel containing its
own haplotypes, which inflates recruitment and accuracy relative to an unseen
sample. `build_truth.py --leave-one-donor-out <dir>` writes donor-excluded
panels; rerunning against those and scoring with `--leave-one-donor-out` gives
the number that generalises. Report both, labelled.

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
