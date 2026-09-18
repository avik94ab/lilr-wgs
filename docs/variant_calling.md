# Variant calling in `lilr-wgs` — what is established and what is not

Companion to [`method.md`](method.md). Scope: sequence variant calling across the
eleven LILR genes from 1000 Genomes 30× srWGS, paralogue similarity, and how a
minimum depth to call a variant should be expressed here.

Nothing in this document changes the copy-number path.

---

## 1. The trust boundary

`PLAN.md` and `docs/method.md` are reliable **through the copy-number
computation**. Downstream of that they describe code that exists but has not been
verified, and the Phase checkboxes should be read as "written," not "validated."

The test suite draws the line precisely:

| module | lines | tests | status |
|---|---|---|---|
| `loci.py` | 343 | `test_loci.py` (14) | verified |
| `loci_chm13.py` | 264 | `test_loci_chm13.py` (39) | verified |
| `coverage.py` | 501 | `test_coverage.py` (27) | verified |
| `depth_model.py` | 312 | `test_depth_model.py` (38) | verified |
| `cn.py` | 625 | `test_cn.py` (37) | verified |
| `realign.py` | 495 | `test_realign.py` (44) | verified |
| `assign.py` | 259 | `test_assign.py` (16) | verified |
| `portable/lilr_cn.py` | — | `test_portable.py` (35) | verified |
| **`callability.py`** | **205** | **none** | **unverified** |
| **`genotype.py`** | **459** | **none** | **unverified** |
| **`sequences.py`** | **359** | **none** | **unverified** |

**1,023 lines carrying the callable track, variant calling and allele naming have
no test coverage**, while everything up to and including assignment has 250
behind it. (`tests/` greps that appear to touch the three modules are incidental:
a string inside a docstring assertion, and a module-name check in the
panel-independence leak test.)

The code is real, not stubbed — `genotype.py` genuinely invokes `HaplotypeCaller`
with `-ploidy` and `-L callable.bed`, `AddOrReplaceReadGroups`,
`whatshap polyphase` and `bcftools consensus`. The problem is not absence; it is
that nothing has checked whether it is right.

**What the 101-donor validation scored was copy number.** LILRB3 100%, LILRA6
99%, LILRA3 95% are copy-number figures. `PLAN.md` §10 says so itself: allele
sequence has not been scored. There is therefore no evidence, at present, for any
claim about variant calls or allele assignments from this pipeline.

### 1.1 One consequence worth isolating

`depth_model.py` is verified and is pure arithmetic on numbers. `callability.py`,
which feeds it real evidence from a BAM, is not. So the thresholds are sound, and
whether they are being applied to correct per-position depth counts is unknown.

`gather_evidence()` is where that risk sits. It walks a pileup with
`stepper="nofilter"`, skips secondary and supplementary alignments, counts
MAPQ ≥ floor separately, and detects shared reads by either the `ZS` tag or the
`shared_names` table. Each of those is a decision that could be individually
wrong without the output looking wrong: a pileup depth that silently differs from
`samtools depth` in deletion/refskip handling, a `min_base_quality=13` filter that
removes more than intended at low coverage, or a shared-read count that misses
reads whose tag did not survive the FASTQ hop and whose name is absent from the
table. Any of those produces a plausible-looking depth and a wrong callable track.

---

## 2. What must happen before any downstream claim

In this order. Each step is cheap relative to the cost of discovering the problem
later.

**2.1 Unit-test the three modules, pure logic first.** Following the repo's
existing separation, `classify()` and `mask_sequence()` in `callability.py` and
the naming logic in `sequences.py` are functions of numbers and strings and can be
tested without a BAM, exactly as `depth_model` and `arbitrate` are. Do this before
touching anything with a cluster in it.

**2.2 Test `gather_evidence()` against an independent depth measurement.** Build a
small fixture BAM with known per-position depth and assert agreement with
`samtools depth -a` position by position, including at deletions, refskips and
secondary alignments. This is the single most valuable test in the set, because it
is the junction between verified arithmetic and unverified I/O.

**2.3 Pin λ₁'s unit conversion with a test.** λ₁ is measured on the CRAM slice and
applied to a realigned per-gene BAM, after extraction, recruitment, arbitration
and realignment have each dropped reads; `recruitment_efficiency()` converts
between the two (0.813 on HG00096).

**The conversion is already applied** — `genotype.py` does it explicitly:

```python
efficiency = float(cov.get("efficiency", 1.0))
lambda1 = float(cov["lambda1"]) * efficiency
```

and the local `gc_lookup` closure scales that already-converted value by the GC
factor. So this is not an outstanding bug. What is outstanding is that **nothing
tests it**, and the failure would be invisible: dropping the efficiency factor
moves every threshold by ~19% and changes no field in any output.

Note also that `genotype.py` reconstructs the conversion inline from the coverage
JSON rather than calling `CoverageModel.lambda_at()`, because it reads a JSON file
rather than holding the object. That is reasonable, but it is a second copy of the
same arithmetic, and CLAUDE.md warns in the other direction ("callers should use
`lambda_at()` rather than `.lambda1` directly"). A test should assert the two
agree, so the copies cannot drift.

**2.4 Then score allele sequence on the 101-donor overlap.**
`validation/build_truth.py` already derives allele sequences from the panels;
`validation/compare_cn.py` is the scoring template. Two requirements carried over
from the copy-number work:

- **Score as-is and leave-one-donor-out.** Copy number cannot be inflated by panel
  circularity because it is measured before a panel is opened — verified byte for
  byte (`PLAN.md` §10). Recruitment *is* inflated, by up to 11% at LILRB3, and
  allele sequence runs through recruitment. The caveat that does not bind copy
  number does bind this.
- **Report per class, never as one aggregate.** `PLAN.md` §14 is explicit that both
  failed LILRA3 windows "would have passed an integer-agreement check on most
  samples," and only the per-class estimate breakdown separated them. The same
  reasoning applies to sequence concordance: report per gene, and within gene by
  region type (gene-unique, shared-block, intronic).

---

## 3. "DP ≥ 6 to confidently call a variant" — the framing is wrong here

This sits in the verified zone, so the analysis below is on firm ground. Both
tables were recomputed from `depth_model.thresholds_for()` at λ₁ = 18.0, which is
the value this project measures (median 17.9–18.6 across the `PLAN.md` §13–14
runs).

`depth_model.py` exists to remove fixed depth constants. The predecessor filters
at `FMT/DP >= max(20, 10*CN)`, and at 30× that "does not remove the doubtful tail,
it removes half the data, and the half it removes is not random — it is the
GC-extreme, the paralogue-adjacent and the repeat-flanked positions, which is to
say exactly the positions where LILR genotypes differ." A flat `DP ≥ 6` is the
same class of object, only smaller.

### 3.1 It is looser than what the model already enforces

| CN | dispersion | expected | floor | ceiling | floor set by |
|---|---|---|---|---|---|
| 1 | Poisson | 18.0 | **8** | 28 | distribution |
| 2 | Poisson | 36.0 | **21** | 51 | distribution |
| 3 | Poisson | 54.0 | **36** | 72 | distribution |
| 1 | r = 10 | 18.0 | **5** | 36 | reads_per_copy |
| 2 | r = 10 | 36.0 | **11** | 61 | distribution |
| 3 | r = 10 | 54.0 | **23** | 85 | distribution |

At every copy number above 1 the existing floor already exceeds 6 by a wide
margin. Adding `DP ≥ 6` as a filter would **admit positions the model currently
rejects**.

### 3.2 Six reads total and six reads per copy are different claims

At CN 2, six reads total gives a true heterozygote a 0.5⁶ ≈ **1.6%** chance of
contributing no alternate read, and **10.9%** of contributing zero or one. Six
reads per haploid copy — twelve at CN 2 — gives 0.5¹² ≈ **0.024%**. The gap widens
with ploidy.

The per-copy form of the constant already exists:

```python
MIN_READS_PER_COPY = 5                      # depth_model.py
per_copy_floor = MIN_READS_PER_COPY * copies
floor = max(distributional_floor, per_copy_floor)
```

### 3.3 The correct expression, and how to justify it

If the intent is six reads per haploid copy, the change is one line —
`MIN_READS_PER_COPY = 6`. Measured effect:

| case | at 5 | at 6 | binds? |
|---|---|---|---|
| CN 1, Poisson | 8 | 8 | no — distribution higher |
| CN 2, Poisson | 21 | 21 | no |
| CN 1, r = 20 | 5 | **6** | yes |
| CN 1, r = 10 | 5 | **6** | yes |
| CN 2, r = 10 | 11 | **12** | yes |
| CN 3–4, any | unchanged | unchanged | no |

A modest tightening at low copy number under real overdispersion, and a no-op
elsewhere — the right shape for a floor whose purpose is to stop a haplotype being
called off a coin flip, and consistent with the existing comment ("three reads is
a coin flip and five is the point where an allele balance starts to mean
something").

Per the repo's convention — record the measurement next to the constant it
justifies — **do not make this change on argument alone.** It cannot be justified
until §2.4 exists: score sequence concordance at both values and put the numbers
in the comment. If they are identical, leave it at 5 and record that.

### 3.4 One threshold, not two

The predecessor masked the consensus at 15 and filtered variants at 20, so a
position with depth in [15, 20) survived the mask but lost its call. Here the
callable track does both jobs and reaches `bcftools` as `-T callable.bed`. Any new
depth requirement enters through `depth_model`, never as a second filter
expression downstream.

---

## 4. Paralogue similarity

**Assignment is verified** (`test_assign.py`, 16 tests) and handles the hard case
deliberately: `assign.arbitrate()` routes each pair to its best gene by paired
alignment score, but a tie confined to one `shared_group` is kept in **every** tied
gene rather than dropped from all. LILRA6 and LILRB3 share a ~4.6 kb block with no
gene-diagnostic 31-mers, and a competitive discard removes it from both at once —
26% of the LILRB3 CDS in the capture cohort. PING does the same for KIR2DL5A/B.
This is the opposite of a margin-based discard and must not be "fixed."

The shared verdict travels in `shared_pairs.tsv` because reads are recruited
against a pangenome panel and called against a single reference, so a FASTQ sits
between them and BAM tags do not survive it.

**The two-sided depth gate is sound in principle and unverified in application.**
`depth_model` rejecting above the upper quantile is the right treatment — in a
cluster of ~90%-identical paralogues a position at three times expected depth is
an unseparated pile-up, and a heterozygote there is a paralogous sequence variant
— and `effective_copies()` correctly stops that backfiring on shared blocks, where
judging against the gene's own copy number would reject ~90% of the LILRA6 CDS.

But whether `gather_evidence()` counts `n_shared` correctly across the FASTQ hop is
exactly the untested question in §1.1, and the failure is directional: if it
**undercounts**, `effective_copies()` receives too low a `shared_fraction`, the
expectation is too low, and the ceiling rejects real shared-block positions as
pile-ups. Test that before trusting a callable fraction at LILRA6 or LILRB3.

**Statuses stay distinct.** `ok` / `low_depth` / `high_depth` / `low_mapq` /
`paralog_ambiguous` / `no_model`, and `measured` / `not_measured` / `failed`. Never
a boolean or a bare `N`: at LILRA3 zero is a common true state, so a failed query
reported as zero manufactures deletions, and a non-ALT-aware GRCh38 alignment
reads near zero at every MAPQ-20 window for everyone alike, which looks exactly
like a cohort of LILRA6 deletion homozygotes.

---

## 5. Open questions that bound the answer

**How much of LILRA6/LILRB3 is callable at 30× at all?** Paralogue-unique sequence
is 2,900 bp and 1,881 bp; outside it, assignment keeps shared-block reads in both
genes. `summarise()` already emits `callable_fraction` per (sample, gene) — the
measurement exists and has not been reported. Report it split by gene-unique,
shared-block and intronic, since averaging hides the structure that matters. This
number bounds what can be claimed about those two genes regardless of how well
anything else works.

**Does 30× support CN > 3 at LILRA6?** At λ₁ ≈ 18, separating CN 4 from CN 5 means
separating 72× from 90×. The CHM13 run reproduced the tails at 173 samples (CN 4:
6, CN 5: 1), which is thin. This is a copy-number question, but it bounds the
ploidy the caller is ever handed.

---

## 6. On the `lilr-atlas` `_full` references

**Not verifiable from this repository** — `lilr-atlas` is not present on the
machine these notes were written on, so everything in this section is carried from
the source document and should be re-checked against the files before being acted
on.

`lilr-atlas/sequences/*_gDNA_full.fasta` is reported to span first to last exon
plus introns plus ±1500 bp, with `_full` minus `_CDS` being exactly 3000 bp for ten
of eleven genes, and each file an allele catalogue of 99–424 records with
three-field names (`LILRA3*001:01:01`) rather than a single reference.

**Check whether they are already here** before adding a reference set.
`resources/gdna/{gene}.fasta` is the 465-sequence HPRC panel used for recruitment,
and `resources/bundle/references/{gene}_named.fa` is the single reference used for
calling. That split is deliberate — a divergent allele attracts its own reads at
recruitment, and variants need one coordinate system at calling — and a third set
would need a reason.

**LILRA6 reportedly breaks the 3000 bp pattern** at 2586 bp. Given it is 47%
ambiguous by the 100-mer criterion and already the hardest locus, resolve whether
the flank is clipped by the duplication boundary or the first record is truncated
before relying on the file.

**The third field is unreachable until §2 is done.** The atlas carries gDNA-level
resolution that targeted short reads could not assign, and WGS spanning introns
could in principle reach it — but `sequences.py`, which does the naming, is among
the untested modules, and kir-mapper's documentation warns that its equivalent
intron mode produces apparent novel alleles that are artifacts of repetitive
regions. A three-field claim requires: `sequences.py` tested, exonic and intronic
concordance scored separately on the overlap set, and every distinguishing
position `ok` in the callable track.

---

## 7. Validation set sizing

The repo already answers the question of how many samples to pilot on. **101 donors
have both an HPRC assembly and a 1000 Genomes CRAM**, and `validation/reports/`
holds `pilot5_cn.txt`, `overlap101_asis.txt` and `overlap101_lodo.txt` (plus
counted variants and the JoGo cross-check). The tiers are 5, 101 and 2,504; use
them rather than constructing a new set.

Two limits of that set:

- **Ploidy coverage at the tails.** Check the copy-number distribution across the
  101 before scoring sequence, and if CN 0 at LILRA3 or CN ≥ 4 at LILRA6 is thin,
  say so rather than reporting an aggregate that never exercised those branches.
- **No power for population-level checks.** The Hardy-Weinberg check in `PLAN.md`
  §11 runs on the full 2,504 for a reason. Excess heterozygosity is the most
  sensitive available detector of unseparated paralogues, and at n = 101 it has
  little power. Do not calibrate a collapse filter on the overlap set.

---

## 8. Recommendations retracted, recorded so they are not reintroduced

- A flat `DP ≥ 6` variant filter — rejected; §3. Belongs in `MIN_READS_PER_COPY`,
  and cannot be justified until sequence scoring exists.
- Dropping alignment-score-ambiguous read pairs on a margin — rejected;
  `arbitrate()` keeps shared-group ties in both genes deliberately (§4).
- Masking uncallable positions to a bare `N` — rejected; the callability
  vocabulary is load-bearing (§4).
- A new `resources/gdna_full/` set and a parallel `wgs_*` script tree —
  unnecessary; `src/lilrwgs/` + `scripts/` + `workflow/` and the existing
  two-reference split already cover it.
- A cohort-wide diagnostics rule consumed downstream — forbidden; the DAG has no
  barrier and `cohort_scale` is a leaf for that reason.
- Treating Phase 6 as complete because `PLAN.md` marks it `[x]` — corrected in §1.
  The three downstream modules have no tests, and the 101-donor scores are
  copy-number scores.
