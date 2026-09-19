"""Turning measured coverage into per-position call thresholds.

This module is the reason the project exists, so it is worth being explicit about
what is wrong with the alternative.

`lilr-genotyper` filters variants at ``FMT/DP >= max(20, 10*CN)`` and masks the
consensus below the same number. Those constants are reasonable for targeted
capture at several hundred ×, where 20 reads is a small fraction of what is
available. At 30× WGS a diploid locus carries about 30 reads, so a threshold of
20 sits near the *middle* of the depth distribution: it does not remove the
doubtful tail, it removes half the data, and the half it removes is not random —
it is the GC-extreme, the paralogue-adjacent and the repeat-flanked positions,
which is to say exactly the positions where LILR genotypes differ.

The replacement is a model with three properties the constants lack.

**It is per-sample.** λ₁, the depth one haploid copy yields, is measured from
control loci in the same CRAM. A library at 25× and one at 38× get different
thresholds, which is the point; a fixed constant silently calls the first one
conservatively and the second one loosely.

**It scales with copy number correctly.** A locus at copy number *k* is modelled
as the sum of *k* independent haploid contributions. Negative binomials are
additive in exactly this way, so μ and the dispersion both scale by *k* and the
threshold follows from the distribution rather than from multiplying a constant.

**It is two-sided.** This is the part with no counterpart in the predecessor. In
a cluster of ~90%-identical paralogues, a position at three times its expected
depth is not better supported — it is a pile-up of reads from a paralogue that
read assignment failed to separate, and a "heterozygous" call there is a
paralogous sequence variant in disguise. A one-sided floor accepts every one of
them. Rejecting above the upper quantile costs a little real data at genuinely
high-coverage positions and buys immunity from a class of error that is
otherwise invisible and systematic.

Everything here is a pure function of numbers, so it is unit-testable without a
BAM, a cluster, or a network.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum

# Below this many reads per haploid copy, a haplotype cannot be distinguished
# from its neighbour whatever the total depth says: three reads is a coin flip
# and five is the point where an allele balance starts to mean something. This
# floor binds at low copy number and low coverage, where the distributional
# bound alone would let through positions that are statistically unremarkable
# and still uncallable.
MIN_READS_PER_COPY = 5

# Two-sided tail mass. 0.005 per side keeps 99% of honest positions and puts the
# CN-2 ceiling at roughly 1.7x expected, which is well below the 2x that an
# unseparated paralogue produces and well above ordinary coverage variation.
DEFAULT_ALPHA = 0.005

# PING's two-stage depth thresholds, adopted here at this project's own values.
#
# PING (Hollenbach lab) runs `setup.minDP <- 8` while building candidate
# genotypes and `final.minDP <- 20` when finalising them, applied flat as
# `vcfDT[DP >= minDP]` — permissive during discovery so a variant is not lost
# before it can be assessed, strict at output. `lilr-genotyper` uses 20 as well.
#
# The values here are 6 and 10, deliberately lower than both, because those
# pipelines call KIR from targeted capture at several hundred x while this one
# calls LILR at 30x, where the same constant is a much larger fraction of the
# available evidence. The distributional floor at this project's measured λ₁
# lands at 13 for a diploid locus (HG00119: λ₁ 18.39, efficiency 0.802,
# dispersion 33.3), so 10 admits the band between them.
#
# **These are absolute floors and they replace the distributional bound when
# supplied, rather than being taken alongside it.** That is a real departure
# from the rest of this module, which exists to remove fixed depth constants,
# and it is worth being explicit rather than quiet about: a flat floor does not
# know the sample's coverage, so on a thin library it admits more than it should
# and on a deep one it is nearly inert. The two-sided ceiling and the
# copy-number scaling of `effective_copies` are untouched, so the paralogue
# pile-up protection that a flat threshold has no opinion about still applies.
MIN_DP_SETUP = 6         # candidate discovery: what HaplotypeCaller may consider
MIN_DP_FINAL = 10        # output: what survives into the VCF and the consensus

# Minimum fraction of reads supporting an allele for it to be called, PING's
# `hetRatio`, at PING's value. This is the filter a depth threshold cannot
# replace: at adequate depth, one or two reads from a 97%-identical paralogue
# produce a false heterozygote that passes any floor, and only an allele-balance
# requirement rejects it. In a cluster where LILRA6 and LILRB3 differ by ~3%,
# that is the more common error than thin coverage.
HET_RATIO = 0.25

# Dispersion is estimated from control loci, where per-base depth varies for real
# reasons (GC, mappability) on top of sampling noise. An estimate from too few
# bases, or one that comes out under-dispersed, falls back to Poisson — which is
# the *narrower* distribution, so the fallback is conservative in the direction
# that matters: it will not widen the acceptance interval on bad evidence.
MIN_BASES_FOR_DISPERSION = 5_000


class Callability(str, Enum):
    """Why a position is or is not usable. Never a bare N.

    The predecessor masked low-depth positions to N and left no record of which
    positions those were, so an N in its output could equally mean "no reads",
    "reads but ambiguous" or "a real deletion". Distinguishing them is most of
    the diagnostic value of the mask.
    """

    OK = "ok"
    LOW_DEPTH = "low_depth"
    HIGH_DEPTH = "high_depth"
    LOW_MAPQ = "low_mapq"
    PARALOG_AMBIGUOUS = "paralog_ambiguous"
    NO_MODEL = "no_model"


@dataclass(frozen=True)
class DepthThresholds:
    """The acceptance interval at one copy number, and what produced it."""

    copies: int
    expected: float
    floor: int
    ceiling: int
    alpha: float
    dispersion: float
    floor_source: str       # "distribution", "reads_per_copy", "absolute"

    def classify(self, depth: float) -> Callability:
        if depth < self.floor:
            return Callability.LOW_DEPTH
        if depth > self.ceiling:
            return Callability.HIGH_DEPTH
        return Callability.OK

    def as_row(self) -> dict:
        return {
            "copies": self.copies,
            "expected": round(self.expected, 2),
            "floor": self.floor,
            "ceiling": self.ceiling,
            "floor_source": self.floor_source,
            "dispersion": round(self.dispersion, 2) if math.isfinite(self.dispersion) else "inf",
        }


def estimate_dispersion(mean: float, variance: float, n_bases: int) -> float:
    """Negative binomial ``r`` from the mean and variance of per-base depth.

    For NB(μ, r), variance is μ + μ²/r, so r = μ²/(v − μ). Returns ``inf``
    (Poisson) when the data are under-dispersed or too sparse to say — the
    narrower distribution, and therefore the conservative fallback.
    """
    if n_bases < MIN_BASES_FOR_DISPERSION or mean <= 0:
        return math.inf
    excess = variance - mean
    if excess <= 0:
        return math.inf
    return (mean * mean) / excess


def _nb_variance(mean: float, dispersion: float) -> float:
    if not math.isfinite(dispersion) or dispersion <= 0:
        return mean            # Poisson
    return mean + (mean * mean) / dispersion


def _normal_quantile(p: float) -> float:
    """Standard normal inverse CDF, Acklam's rational approximation.

    Accurate to ~1e-9 over the whole range, which is far more than needed for a
    tail at 0.005 — and it keeps this module dependency-free, so the thresholds
    can be recomputed anywhere, including inside a plotting notebook or a test
    that has no scipy.
    """
    if not 0.0 < p < 1.0:
        raise ValueError(f"quantile needs 0 < p < 1, got {p}")

    a = (-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00)
    b = (-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01)
    c = (-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00)
    d = (7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00)
    p_low, p_high = 0.02425, 1 - 0.02425

    if p < p_low:
        q = math.sqrt(-2 * math.log(p))
        return (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) / \
               ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1)
    if p > p_high:
        q = math.sqrt(-2 * math.log(1 - p))
        return -(((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) / \
                ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1)
    q = p - 0.5
    r = q * q
    return (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5]) * q / \
           (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1)


def thresholds_for(copies: int, lambda1: float, dispersion: float = math.inf,
                   alpha: float = DEFAULT_ALPHA,
                   min_dp: int | None = None) -> DepthThresholds:
    """The depth interval a position at *copies* copies should fall in.

    Args:
        copies: copy number at this locus for this sample. Zero means the gene is
            absent, and the returned interval rejects everything — a position
            with depth at a deleted locus is evidence of mis-assignment, not of a
            variant.
        lambda1: expected depth from one haploid copy, from the coverage model.
        dispersion: negative binomial ``r`` for one haploid copy. Additive, so a
            locus at *k* copies uses ``k * dispersion``.
        alpha: tail mass excluded on each side.

    Both bounds use a normal approximation to the negative binomial. At the means
    involved here — 15 and up, since anything lower is uncallable regardless —
    the approximation's error is far smaller than the uncertainty in λ₁ itself,
    and the exact quantile would add a scipy dependency for a correction that
    changes a threshold by less than one read.
    """
    if copies <= 0:
        return DepthThresholds(copies=0, expected=0.0, floor=1, ceiling=0,
                               alpha=alpha, dispersion=dispersion,
                               floor_source="absent")

    mean = copies * lambda1
    r_k = dispersion * copies if math.isfinite(dispersion) else math.inf
    sd = math.sqrt(_nb_variance(mean, r_k))
    z = _normal_quantile(1 - alpha)

    distributional_floor = mean - z * sd
    per_copy_floor = MIN_READS_PER_COPY * copies

    if min_dp is not None:
        # PING-style absolute floor: flat, and it *replaces* the modelled bound
        # rather than joining it, so asking for 10 gives 10 rather than the
        # higher of 10 and whatever the distribution wanted. The ceiling below
        # is left alone, so excess depth is still rejected.
        floor_value = float(min_dp)
        floor_source = "absolute"
    else:
        floor_value = max(distributional_floor, per_copy_floor)
        floor_source = ("reads_per_copy" if per_copy_floor >= distributional_floor
                        else "distribution")

    return DepthThresholds(
        copies=copies,
        expected=mean,
        floor=max(1, math.ceil(floor_value)),
        ceiling=math.floor(mean + z * sd),
        alpha=alpha,
        dispersion=dispersion,
        floor_source=floor_source,
    )


@dataclass(frozen=True)
class PositionCall:
    """One position's verdict, with the numbers that produced it."""

    depth: float
    effective_copies: float
    status: Callability
    expected: float
    floor: int
    ceiling: int

    @property
    def usable(self) -> bool:
        return self.status is Callability.OK


def effective_copies(copies: int, shared_fraction: float,
                     paralog_copies: int) -> float:
    """Copies' worth of reads expected at a position, allowing for shared blocks.

    Where two paralogues cannot be separated by short reads — LILRA6 with
    LILRB3, LILRB1 with LILRB4 — read assignment deliberately keeps a tied pair
    in *both* genes rather than discarding it from both, because a competitive
    discard deletes the shared block from both at once. The consequence is that
    depth in such a block counts reads from both genes, so the expectation there
    is the pair's combined copy number, not the gene's own.

    Without this correction the two-sided gate would reject the entire shared
    block as over-covered — which is to say it would reject 90% of the LILRA6
    coding sequence in every sample, and do it for a reason that is a modelling
    error rather than a property of the data.

    Args:
        shared_fraction: fraction of reads at the position carrying the shared
            tag, from 0 (gene-unique) to 1 (entirely shared).
    """
    shared_fraction = min(max(shared_fraction, 0.0), 1.0)
    return copies + shared_fraction * paralog_copies


def call_position(depth: float, copies: int, lambda1: float,
                  *, dispersion: float = math.inf, alpha: float = DEFAULT_ALPHA,
                  shared_fraction: float = 0.0, paralog_copies: int = 0,
                  mapq_fraction: float = 1.0, min_mapq_fraction: float = 0.5,
                  min_dp: int | None = None) -> PositionCall:
    """Classify one position.

    Order matters. MAPQ is checked first because a position whose reads are
    nearly all multi-mapping is uninterpretable regardless of how many there are,
    and reporting it as a depth problem would send someone looking at coverage
    when the issue is that the sequence is not unique.

    Args:
        mapq_fraction: fraction of reads at the position above the MAPQ floor.
        min_mapq_fraction: below this, the position is LOW_MAPQ. 0.5 rather than
            something stricter because in this cluster a position where half the
            reads place confidently still carries usable signal, and the
            shared-block machinery above handles the rest.
    """
    if lambda1 <= 0:
        return PositionCall(depth, 0.0, Callability.NO_MODEL, 0.0, 0, 0)

    # `depth > 0` guards the MAPQ test, and it is not a formality. At zero depth
    # `mapq_fraction` is 0/0, which PositionEvidence reports as 0.0 — so without
    # this an uncovered position fails the MAPQ test and is labelled LOW_MAPQ
    # despite there being no reads to have a mapping quality. The whole point of
    # separating these statuses is that they prescribe different actions: a gene
    # that is N for low depth needs more coverage, one that is N for low MAPQ
    # never will. A gene at copy number 0 — a LILRA3 deletion homozygote, a
    # LILRA6 null — has no reads anywhere in it, and used to report its entire
    # length as LOW_MAPQ, which is the most misleading answer available.
    if depth > 0 and mapq_fraction < min_mapq_fraction:
        return PositionCall(depth, 0.0, Callability.LOW_MAPQ, 0.0, 0, 0)

    eff = effective_copies(copies, shared_fraction, paralog_copies)
    # Thresholds are defined on integer copies; interpolating the interval would
    # imply a precision the model does not have, so round and carry the exact
    # effective value in the result for anyone who wants it.
    t = thresholds_for(max(1, round(eff)), lambda1, dispersion, alpha,
                       min_dp=min_dp)
    status = t.classify(depth)

    # A shared block that is over-covered even after accounting for the paralogue
    # is a different finding from an ordinary pile-up: it says assignment let in
    # reads from somewhere the model does not know about.
    if status is Callability.HIGH_DEPTH and shared_fraction > 0.5:
        status = Callability.PARALOG_AMBIGUOUS

    return PositionCall(depth, eff, status, t.expected, t.floor, t.ceiling)


def summarise(calls: list[PositionCall]) -> dict:
    """Counts by status, plus the callable fraction.

    The callable fraction is the headline number for a (sample, gene): it says
    how much of the gene the data could actually speak to, which is the question
    an N-heavy consensus leaves unanswered.
    """
    total = len(calls)
    counts = {s.value: 0 for s in Callability}
    for c in calls:
        counts[c.status.value] += 1
    return {
        "n_positions": total,
        "callable_fraction": round(counts[Callability.OK.value] / total, 4) if total else 0.0,
        **{f"n_{k}": v for k, v in counts.items()},
    }
