# Cross-validation against JoGo-LILR

The HPRC overlap ([`README.md`](README.md)) scores this pipeline against
assemblies. This scores it against **another caller**, which is a different kind
of evidence: a shared blind spot between our depth model and an assembly-derived
truth set is conceivable, whereas one shared with a method that starts from
CNVnator bins and ends in a clustering step is not.

[JoGo-LILR](https://doi.org/10.1016/j.humimm.2025.111272) (Nagasaki et al., *Hum
Immunol* 2025;86(3):111272) is the published LILRB3/LILRA6 copy-number caller for
srWGS, distributed as `JoGo-LILR_v1.1.tgz` from
<https://nagasakilab.csml.org/data/JoGo-LILR_v1.1.tgz> (7 GB). It ships two things
that make it usable as a reference rather than only as a method: four-region
normalised coverage for **3,202 1000 Genomes samples**
(`input/LILRB3_LILRA6.cnv.hapmap3202.reference.tsv`) and the paper's own calls for
the same samples (`test/paper.hapmap3202.allele_stable.tsv`). Any 1000 Genomes
cohort we run is already in both.

## What it does differently

Four regions (`input/LILRB3_LILRA6.region.bed`):

| region | GRCh38 |
|---|---|
| LILRB3 core | chr19:54,216,000-54,219,900 |
| LILRB3 gene | chr19:54,216,000-54,223,100 |
| LILRA6 core | chr19:54,236,600-54,239,600 |
| LILRA6 gene | chr19:54,236,000-54,243,000 |

From those it forms `X = LILRB3gene + LILRA6gene` — the pair total — and
`Y = LILRB3core / LILRA6core` — the split, clamped at 4 where LILRA6 is null and
the denominator is 0. It then clusters the cohort by single linkage against the
3,202-sample background and snaps each cluster's centroid to the nearest
theoretical anchor (`CN5_B2A3`: total 5, B=2, A=3).

Two differences from this pipeline are the point of the exercise. It is
**cohort-relative** — a sample is called by where it sits among other samples,
which is the barrier `lilr-wgs` exists to avoid — and it reads the **pair total
and ratio**, not the paralogue-unique window, so it fails differently.

Preprocessing is CNVnator 0.4.1 with `-chrom chr19`, bin 100: **normalisation is
chr19-wide, so our 550 kb slices cannot feed it as-is.** Substituting λ₁ for the
chr19 mean reproduces their region values closely enough to run their caller on
our slices: on `kgp100b`, X correlates at r = 0.994 with slope 1.010 and Y at
r = 0.999, with a systematic offset of +0.167 in X.

## Result

Two cohorts of 100, drawn disjoint from each other and from the 101 HPRC
validation donors (`config/manifest.kgp100.tsv`, seed 20260917;
`config/manifest.kgp100b.tsv`, seed 202609172), against the paper's published
LILRA6 copy number for the same samples:

| route | kgp100 | kgp100b |
|---|---|---|
| **lilr-wgs vs their published calls** | **100/100** | **100/100** |
| their caller on their coverage, vs their published calls | 100/100 | 100/100 |
| their caller on our slices (λ₁-normalised), vs published | 100/100 | 98/100 |

**200/200 on LILRA6**, across CN 0 to 5, in two cohorts drawn a day apart.

The two misses in the last row are the λ₁-substitution route, not this pipeline:
HG00178 (called 5, published 4) and HG02805 (3, published 2), both one copy high,
both consistent with the +0.167 offset in X pushing a borderline sample over an
anchor boundary. `lilr-wgs` called both correctly and flagged both as low
confidence (0.48 and 0.33), which is what measuring the unique window directly
rather than going through the pair total buys.

Outputs: `results/kgp100/lilra6_cn_jogo.tsv` and
`results/kgp100b/lilra6_cn_jogo.tsv` (under `results/`, so untracked). Both carry
per-sample X and Y from each route beside the calls.

## What this does not establish

**LILRB3 agreement is much weaker evidence than the LILRA6 figure**, and the
100/100 there should not be quoted alongside it. Their `diploid_stable` anchor
grid contains only B1/B2 haplotype types, so LILRB3 CN 2 is very nearly what the
method is able to return; agreeing with it is close to agreeing that most people
have two copies of LILRB3, which both methods would manage without being right
about anything.

The comparison is also on **copy number only**. Neither the paper's calls nor
ours are haplotype sequences, and the allele-type strings the two methods emit
are not in the same vocabulary.

Finally, the reference table is 1000 Genomes. Both cohorts are drawn from it, so
this says nothing about a cohort whose ancestry composition is different from
theirs — which for a cohort-relative method is a live question, and for this one
is the thing λ₁ is supposed to make moot.
