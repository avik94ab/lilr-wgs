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

Every sample in the 1000 Genomes 30x collection — **2,504**, all of them in their
3,202-sample published table — against the paper's LILRA6 copy number:

**2,497/2,504 = 99.7% agreement**, with no sample refused by either side.

Split by whether this pipeline flagged the call as ambiguous, which is what the
flag is for:

| | n | agree |
|---|---|---|
| confident | 2,125 | 2,123 (99.9%) |
| flagged | 379 | 374 (98.7%) |

All seven disagreements sit at CN >= 3, where 30x separates adjacent copy numbers
by ~15x of depth, and five of the seven are flagged:

| sample | theirs | ours | our estimate | confidence | flagged |
|---|---|---|---|---|---|
| HG00238 | 3 | 4 | 3.602 | 0.203 | yes |
| HG00315 | 3 | 4 | 3.508 | 0.017 | yes |
| HG00742 | 7 | 6 | 6.975 | 0.000 | yes |
| HG01187 | 3 | 4 | 3.645 | 0.290 | yes |
| HG01619 | 5 | 4 | 4.084 | 0.832 | no |
| NA18933 | 5 | 6 | 5.864 | 0.727 | no |
| NA20827 | 5 | 6 | 5.539 | 0.079 | yes |

NA20827 is the one sample that also disagreed with the HPRC assembly truth set,
in the same direction and by the same amount: truth 5, ours 6 (`PLAN.md` §9).
Two methods that share no failure mode now say 5, so that is a genuine miss and
not a truth-set artefact — the only one either exercise has produced.

Before the full cohort was run, the same comparison was made on two disjoint
100-sample draws that share no donor with the 101 HPRC validation donors
(`config/manifest.kgp100.tsv`, seed 20260917; `config/manifest.kgp100b.tsv`, seed
202609172), which is the version to cite where the panels matter:

| route | kgp100 | kgp100b |
|---|---|---|
| **lilr-wgs vs their published calls** | **100/100** | **100/100** |
| their caller on their coverage, vs their published calls | 100/100 | 100/100 |
| their caller on our slices (lambda1-normalised), vs published | 100/100 | 98/100 |

The two misses in the last row are the lambda1-substitution route, not this
pipeline: HG00178 (called 5, published 4) and HG02805 (3, published 2), both one
copy high, both consistent with the +0.167 offset in X pushing a borderline
sample over an anchor boundary. `lilr-wgs` called both correctly and flagged both
as low confidence (0.48 and 0.33), which is what measuring the unique window
directly rather than going through the pair total buys.

Outputs: `validation/reports/jogo_kgp2504.txt` for the full cohort;
`results/kgp100/lilra6_cn_jogo.tsv` and `results/kgp100b/lilra6_cn_jogo.tsv` for
the two 100-sample routes, which carry per-sample X and Y from each route beside
the calls. Reproduce the full-cohort report with:

```bash
python validation/compare_jogo.py results/kgp2504/cn_calls.tsv \
    --published $JOGO/test/paper.hapmap3202.allele_stable.tsv \
    -o validation/reports/jogo_kgp2504.txt
```

## What this does not establish

**LILRB3 agreement is much weaker evidence than the LILRA6 figure** and should
not be quoted alongside it. Their `diploid_stable` anchor grid contains only
B1/B2 haplotype types, so LILRB3 CN 2 is very nearly the only answer the method
can return — and this pipeline calls 2,500 of 2,504 samples CN 2 as well.
Agreeing there is close to agreeing that most people have two copies of LILRB3,
which both methods manage without being right about anything. `compare_jogo.py`
therefore scores LILRA6 and declines to score LILRB3 at all.

The comparison is also on **copy number only**. Neither the paper's calls nor
ours are haplotype sequences, and the allele-type strings the two methods emit
are not in the same vocabulary.

Finally, their background table *is* 1000 Genomes, and the comparison is now the
whole of it. So this says nothing about a cohort whose ancestry composition
differs from theirs — a live question for a cohort-relative method, and the thing
λ₁ is meant to make moot for this one. The honest reading of 99.7% is that the
two methods agree on the samples their method was calibrated on.
