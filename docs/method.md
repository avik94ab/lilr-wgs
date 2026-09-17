# Method

What `lilr-wgs` measures, how, and which claims rest on measurement rather than
assumption. Where a number appears here it was measured on real data and the
sample is named; where something is assumed, it says so.

## 1. The problem with the constants

`lilr-genotyper` filters variants at `FMT/DP >= max(20, 10·CN)` and masks the
consensus below the same threshold. For targeted capture at several hundred ×
those are reasonable: 20 reads is a small fraction of what is available, and the
positions they exclude are genuinely doubtful.

At 30× WGS a diploid locus carries about 30 reads. A threshold of 20 is then not
a filter but a mask across the middle of the distribution — and the positions it
removes are not a random 40%. They are the GC-extreme, repeat-flanked and
paralogue-adjacent ones: precisely where LILR haplotypes differ, and precisely
what the pipeline exists to resolve.

So every threshold here is derived rather than fixed.

## 2. λ₁ — what one haploid copy looks like

λ₁ is the depth a single haploid copy of unique sequence yields in this sample.
It is measured from control loci, not from other LILR genes, and that is the
structural advantage srWGS has over capture. In a capture panel the only
available reference for "one copy" is another gene in the same panel, so copy
number is inherently relative, a cohort must be called as one batch, and a
panel-wide efficiency shift is invisible. Here the controls are ordinary genomic
sequence a few hundred kilobases away in the same CRAM.

Controls, all measured clean by a 100-mer uniqueness criterion:

| locus | interval (GRCh38) | span | ambiguous | inside the alt placement |
|---|---|---|---|---|
| PRPF31 | chr19:54,115,410–54,131,719 | 16.3 kb | 1.7% | yes |
| TMC4 | chr19:54,160,095–54,173,250 | 13.2 kb | 0.0% | yes |
| MBOAT7 | chr19:54,173,412–54,189,882 | 16.5 kb | 0.4% | yes |
| TTYH1 | chr19:54,415,219–54,436,904 | 21.7 kb | 0.1% | yes |
| PRKCG | chr19:53,882,196–53,907,652 | 25.5 kb | 0.4% | no |
| CACNG6 | chr19:53,991,148–54,012,666 | 21.5 kb | 0.0% | no |
| PPP6R1 | chr19:55,229,778–55,259,017 | 29.2 kb | 0.2% | no |

LAIR1 sits closer to the targets and was the obvious fifth, but it is 28.4%
ambiguous against LAIR2 and is not used. The KIR cluster from chr19:54.72 Mb is
excluded outright: it is the most copy-number-variable region in the genome.

λ₁ is half the median control depth. The median across loci, not the mean — one
control overlapping an unannotated CNV in one sample should move the estimate by
nothing, and with four loci a mean would move it by a quarter of the error.

Measured on HG00096: control depths 33.9–38.7, λ₁ = 19.08.

### Dispersion

Per-base depth at a control varies for real reasons — GC, mappability — on top of
sampling noise, so it is overdispersed relative to Poisson. Fitting NB(μ, r) from
the moments gives r, and because a diploid locus is two haploid copies and NB is
additive, the per-haploid r is half the fitted value. HG00096: r = 22.0.

An under-dispersed or poorly-supported estimate falls back to Poisson, which is
the *narrower* distribution — the fallback is conservative in the direction that
matters, since it will not widen the acceptance interval on weak evidence.

### GC

The LILR genes are not at the genomic mean GC, so an uncorrected λ₁ is
systematically wrong at exactly the loci of interest, in a gene-dependent
direction. The curve is fitted from the control loci: each GC bin's own level is
a **median** (so one pile-up does not drag a bin) and the normaliser is the
**mean** over all control bases, because λ₁ is itself a mean and the two have to
compose. Normalising by the median instead leaves a silent, GC-dependent bias in
every threshold — a unit test enforces the composition.

## 3. Two units for λ₁

λ₁ is measured on the CRAM slice. It is *applied* to a per-gene BAM produced much
later: extract to FASTQ → recruit against a pangenome panel → arbitrate
cross-mapping → realign to a single per-locus reference. Each step drops reads,
so the two numbers are in different units.

Measured on HG00096: the CRAM says a diploid locus should carry 38.2, and
LILRB1's realigned median was 31 — 19% down. Applied without correction, that
shortfall put 27% of the gene below a floor it should have cleared, which is the
predecessor's failure mode arriving by a different route.

`recruitment_efficiency()` measures the factor from separable, copy-stable anchor
genes (LILRA1, LILRA2, LILRB2, LILRB5) realigned along the same path. HG00096:
anchor medians 27–33, factor 0.813, effective λ₁ = 15.5.

This does not reintroduce a relative baseline. Copy number still comes from the
CRAM-level measurement against external controls; what is estimated here is the
efficiency of a fixed software path, which is a property of the pipeline and not
of the sample's biology. The anchors are copy-stable at 2 in 231 of 232 HPRC
donors, which the truth set confirms rather than assumes.

## 4. Callability

A position at copy number *k* is modelled as NB(k·λ₁, k·r), and is callable when
its depth falls inside a central interval of that distribution.

**The floor** is the greater of the distributional bound at α = 0.005 and
5 reads per haploid copy. The second term binds at low coverage, where the
distributional bound goes to nothing but a haplotype still needs reads on it to
be distinguished from its neighbour.

**The ceiling** is the distributional bound on the other side, and has no
counterpart in the predecessor. In a cluster of ~90%-identical paralogues a
position at three times its expected depth is a pile-up of reads that assignment
failed to separate, and a heterozygous call there is a paralogous sequence
variant wearing a heterozygote's clothes. A one-sided floor accepts every one of
them, silently and systematically.

**Shared blocks are the exception that makes the ceiling workable.** Where
LILRA6/LILRB3 reads are deliberately kept in both genes, depth counts both, so
the expectation is the pair's *combined* copy number:
`effective_copies = own + shared_fraction · paralog`. Judged against the gene's
own copy number, the gate would reject ~90% of the LILRA6 coding sequence in
every sample — for a modelling reason rather than a data one.

The same interval bounds both the variant call (as `-L`/`-T` to GATK and
bcftools) and the consensus mask. In the predecessor these drifted apart — a mask
at 15 under a filter at 20 — so a position with depth in between survived the
mask but could never receive a variant call, and `bcftools consensus` emitted the
reference base: an unflagged reference-biased call rather than an honest N.
Deriving both from one track makes that gap unrepresentable.

Every position's verdict is written to a track with the numbers behind it, and
the reasons are distinct: `low_depth`, `high_depth`, `low_mapq`,
`paralog_ambiguous`, `no_model`. A gene that is 40% N from low depth needs more
coverage; one that is 40% N from paralogue ambiguity never will be.

Measured on HG00096 at λ₁_eff = 15.5:

| gene | callable | low depth | high depth | low MAPQ | paralog | phased |
|---|---|---|---|---|---|---|
| LILRA1 | 95.6% | 307 | 60 | 53 | 0 | yes |
| LILRA4 | 96.9% | 234 | 15 | 18 | 0 | yes |
| LILRB1 | 80.1% | 1735 | 0 | 13 | 0 | yes |
| LILRA6 | 78.5% | 646 | 0 | 761 | 54 | yes |
| LILRB3 | 75.6% | 1491 | 1 | 644 | 32 | yes |

The separable genes are ~96% callable, and the loss at LILRA6/LILRB3 is
concentrated in MAPQ and paralogue ambiguity rather than depth — which is where
the sequence says it should be.

## 5. Copy number

Three genes are called from data; the other eight are fixed at 2.

**LILRA6 and LILRB3** — depth over paralogue-unique windows (2,900 bp and
1,881 bp at the 3' ends, where an activating and an inhibitory receptor genuinely
differ), divided by λ₁. These windows exist only at MAPQ 20, so the measurement
is conditional on §6.

**LILRA3** — not on the GRCh38 primary assembly at all; the reference chromosome
carries the common ~6.7 kb deletion. Read twice:

- *depth* over the four alt contigs carrying it, at MAPQ 0 (those reads multi-map
  four ways by construction, so the MAPQ floor is deliberately dropped), summed
  and divided by **one** interval's length rather than four;
- *the junction* at chr19:54,296,977, where a LILRA3-bearing chromosome's reads
  soft-clip and a deleted one's cross cleanly. `clipped/(clipped+spanning)`
  estimates half the copy number.

The two share no failure mode — depth dies if the CRAM's reference lacks the alt
contigs, the junction is weakened if the aligner had them and moved the clipped
reads there — so they are both reported and their disagreement is a field.
HG00096: depth 0.064 copies, junction 0 of 19 reads clipped. Both say CN 0.

**The pair check.** LILRA6 is also obtainable as pooled LILRA6+LILRB3 depth at
MAPQ 0 minus LILRB3 from its unique window. Reads from both genes multi-map
across both gene bodies, so mean depth over the union is the pair's copy number
divided by two — the span weight matters, and omitting it makes every 2+2 sample
imply zero LILRA6 copies. HG00096: pooled 3.98, implied LILRA6 1.98 against 2.00
from the unique window.

Copy number is **absolute and per-sample**. `refine_cohort()` checks whether the
cohort's estimates cluster on a unit of 1.0 and reports a systematic offset
rather than dividing it out, so nothing downstream depends on a cohort fit.

## 6. Was the alignment ALT-aware?

GRCh38 carries nine LRC alt haplotypes, and UCSC places all nine at one primary
interval — **chr19:54,025,633–55,084,318** — which contains every window used
here, targets and controls alike.

Given the alt contigs in the index *and* the `.alt` file naming them, bwa computes
MAPQ over primary hits alone and a read in an LRC window keeps it. Given the alt
contigs without the `.alt` file, that read has nine more equally good placements
and comes back MAPQ 0. Asking the 100-mer uniqueness question with the alts in
the search space: LILRA6 and LILRB3 go from 39% and 54% ambiguous to **100%**,
and the four LRC controls from 0.0–1.7% to 96–99%.

So in a non-ALT-aware alignment there is nothing left — and, critically, the
result is not an error. Every MAPQ-20 count reads near zero, and divided by a
live baseline that is a confident zero copies. A cohort of those is a cohort of
LILRA6 deletion homozygotes: a wrong answer that arrives looking like a right
one, at a gene where deletion is real and interesting.

Three loci outside the placement make the case distinguishable, and the pipeline
measures it per sample:

- `q20_lrc` — MAPQ-20 over MAPQ-0 depth at the LRC controls;
- `q20_outside` — the same outside the placement, which separates "the alt
  contigs are eating MAPQ" from "MAPQ is degraded everywhere";
- `dilution` — MAPQ-0 depth outside over inside.

When `q20_lrc` falls below 0.30, LILRA6 and LILRB3 are reported as
**`not_measured`**, not as zero.

Measured on HG00096: q20_lrc = 0.996, q20_outside = 0.999, dilution = 0.930 →
`alt_aware`. The CRAM header confirms the mechanism (bwa-mem 0.7.15 against
`GRCh38_full_analysis_set_plus_decoy_hla.fa`, all nine alt contigs in the
header), but the verdict is measured rather than inferred from the header,
because the `.alt` file's presence is not recorded there.

## 7. Duplicates

The predecessor skips MarkDuplicates deliberately: in targeted capture, pairs
legitimately share start coordinates because probe geometry drives where reads
begin, and marking flagged 70–80% of pairs and crushed depth to nothing.

In WGS a shared start coordinate means what it usually means. The 1000 Genomes
CRAMs arrive duplicate-marked by NYGC (picard 2.4.1, visible in the header), so
the flags are honoured at the slice and nothing is recomputed — recomputing them
on a 550 kb slice would be wrong anyway, since duplicate detection needs the whole
library to judge.

## 8. Validation

The HPRC panels the pipeline aligns against carry their donor in the FASTA
header, so they are also a truth set: 232 donors, of whom **101 have a 1000
Genomes 30× CRAM**.

Truth copy-number distribution from the panels:

| gene | distribution |
|---|---|
| LILRA3 | 0: 24, 1: 64, 2: 144 |
| LILRA6 | 0: 1, 1: 13, 2: 141, 3: 54, 4: 16, 5: 5, 6: 2 |
| LILRB3 | 1: 3, 2: 228, 3: 1 |
| LILRA1/A2/A5/B2/B5 | 2: 231, 3: 1 |

The LILRA3 CN 0 entries are inferred: a donor present in at least 8 of the 11
panels is taken to have a well-resolved LRC, so absence from one more means
deletion rather than a failed assembly. That gives a deletion allele frequency of
24%, consistent with published estimates. The inference is recorded per row so it
can be scored with or without.

The stability of LILRA1, LILRA2, LILRB2 and LILRB5 at 2 copies in 231 of 232
donors is what justifies using them as the recruitment-efficiency anchors — a
measurement, not an assumption.

**The circularity.** A donor in the overlap is aligned against a panel containing
its own haplotypes, which inflates recruitment and accuracy relative to an unseen
sample. `build_truth.py --leave-one-donor-out` writes donor-excluded panels, and
that score is the one that generalises. Both are reported, labelled.

## 9. Pilot result

Five donors from the overlap, chosen to span the hard range rather than to be
easy — LILRA6 at CN 1, 2, 2, 3 and 4, and LILRA3 at CN 0 and 2:

| sample | gene | truth | called | estimate |
|---|---|---|---|---|
| HG02155 | LILRA6 | 1 | 1 | 0.961 |
| HG00099 | LILRA6 | 2 | 2 | 2.336 |
| NA18608 | LILRA6 | 2 | 2 | 2.093 |
| HG00140 | LILRA6 | 3 | 3 | 3.149 |
| HG02922 | LILRA6 | 4 | 4 | 3.779 |
| NA18608 | LILRA3 | 0 | 0 | 0.002 |
| HG00099 / HG00140 / HG02155 / HG02922 | LILRA3 | 2 | 2 | 2.047 / 1.933 / 2.034 / 1.989 |
| all five | LILRB3 | 2 | 2 | 1.758 – 2.267 |

**15 of 15 correct.** All five samples read `alt_aware`; λ₁ ranged 16.6–20.6 and
the recruitment efficiency 0.750–0.799.

This is a pilot, and it should be read as one. Five samples cannot distinguish a
method that is right from one that is right on easy cases, the scoring is *as-is*
rather than leave-one-donor-out, and LILRA6 at CN 4 is represented once. What it
does establish is that the measurement is not grossly biased and that the hard
classes — a LILRA6 hemizygote, a LILRA6 at four copies, a LILRA3 deletion
homozygote — are separable at 30× at all.

It is also the step that found the LILRA3 bug. Before the truth set existed, a
uniform "LILRA3 CN 0" across every sample looked entirely plausible.

## 10. What is not yet established

- Accuracy across the full 101-donor overlap, as-is and leave-one-donor-out.
- Whether 30× supports LILRA6 CN ≥ 4. Separating CN 4 from CN 5 means separating
  60× from 75×: feasible per position, marginal per sample. The truth set
  contains 16 donors at CN 4, 5 at CN 5 and 2 at CN 6.
- Whether the shared-block machinery yields usable *haplotype sequence* at
  LILRA6/LILRB3, as opposed to usable depth. 78.5% and 75.6% callable on one
  sample is encouraging and is not the same claim.
- Whether the negative binomial is the right family for real LILR depth. The unit
  tests check the arithmetic against simulation from the same distribution, which
  says nothing about the fit.
