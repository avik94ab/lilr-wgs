"""Tests for the depth model.

The properties checked here are the ones that would make the model wrong in a
way that produced plausible-looking output: thresholds that do not move with
coverage, an interval that does not scale with copy number, a shared block that
gets rejected as over-covered, a floor that reproduces the predecessor's
behaviour of masking most of a 30x gene.
"""

from __future__ import annotations

import math

import pytest

from lilrwgs.depth_model import (
    Callability,
    DEFAULT_ALPHA,
    MIN_READS_PER_COPY,
    call_position,
    effective_copies,
    estimate_dispersion,
    summarise,
    thresholds_for,
    _normal_quantile,
)

# A 30x library: one haploid copy yields ~15 reads.
LAMBDA1 = 15.0
# Control loci in real data are overdispersed; r ~ 10 is a plausible value and
# the tests that care about the exact number say so.
DISPERSION = 10.0


class TestNormalQuantile:
    @pytest.mark.parametrize("p,expected", [
        (0.5, 0.0),
        (0.975, 1.959964),
        (0.995, 2.575829),
        (0.005, -2.575829),
        (0.001, -3.090232),
    ])
    def test_known_values(self, p, expected):
        assert _normal_quantile(p) == pytest.approx(expected, abs=1e-5)

    def test_rejects_out_of_range(self):
        for bad in (0.0, 1.0, -0.1, 1.5):
            with pytest.raises(ValueError):
                _normal_quantile(bad)


class TestDispersion:
    def test_poisson_when_variance_equals_mean(self):
        assert estimate_dispersion(mean=30, variance=30, n_bases=100_000) == math.inf

    def test_underdispersed_falls_back_to_poisson(self):
        # Narrower than Poisson should not widen the interval on the strength of
        # an estimate that cannot be right.
        assert estimate_dispersion(mean=30, variance=12, n_bases=100_000) == math.inf

    def test_too_few_bases_falls_back_to_poisson(self):
        assert estimate_dispersion(mean=30, variance=90, n_bases=100) == math.inf

    def test_recovers_r_from_moments(self):
        # NB(mu=30, r=10) has variance 30 + 900/10 = 120.
        assert estimate_dispersion(30, 120, 100_000) == pytest.approx(10.0)


class TestThresholds:
    def test_interval_brackets_the_expectation(self):
        t = thresholds_for(2, LAMBDA1, DISPERSION)
        assert t.floor < t.expected < t.ceiling
        assert t.expected == pytest.approx(30.0)

    def test_scales_with_copy_number(self):
        two = thresholds_for(2, LAMBDA1, DISPERSION)
        four = thresholds_for(4, LAMBDA1, DISPERSION)
        assert four.expected == pytest.approx(2 * two.expected)
        assert four.floor > two.floor
        assert four.ceiling > two.ceiling

    def test_scales_with_coverage(self):
        """The whole point of a per-sample model: a deeper library gets a higher
        floor, and a shallower one is not held to a threshold it cannot meet."""
        shallow = thresholds_for(2, 10.0, DISPERSION)
        deep = thresholds_for(2, 25.0, DISPERSION)
        assert deep.floor > shallow.floor
        assert deep.ceiling > shallow.ceiling

    def test_floor_is_far_below_the_predecessors_constant(self):
        """The regression this project exists to prevent.

        `lilr-genotyper` used DP>=20 at CN 2. At 30x that is two thirds of the
        expected depth and masks most of the gene. Whatever this model does, it
        must not land there.
        """
        t = thresholds_for(2, LAMBDA1, DISPERSION)
        assert t.floor < 20, f"floor {t.floor} reproduces the capture-era constant"

    def test_per_copy_floor_binds_at_low_coverage(self):
        """At low λ₁ the distributional bound goes to nothing, and a haplotype
        still needs reads on it to be called."""
        t = thresholds_for(2, 4.0, DISPERSION)
        assert t.floor >= MIN_READS_PER_COPY * 2
        assert t.floor_source == "reads_per_copy"

    def test_distribution_binds_at_high_coverage(self):
        t = thresholds_for(2, 60.0, DISPERSION)
        assert t.floor_source == "distribution"

    def test_poisson_interval_is_narrower_than_overdispersed(self):
        poisson = thresholds_for(2, LAMBDA1, math.inf)
        overdisp = thresholds_for(2, LAMBDA1, DISPERSION)
        assert poisson.ceiling < overdisp.ceiling
        assert poisson.floor > overdisp.floor

    def test_absent_locus_has_nothing_callable(self):
        """CN 0 is a common true answer at LILRA3, not an error state.

        Nothing is callable there, but the two ways of being uncallable say
        different things and must not be merged: no reads is the gene being
        absent as expected, while reads at a locus called absent is a
        contradiction between the copy number and the alignment.
        """
        t = thresholds_for(0, LAMBDA1, DISPERSION)
        assert t.floor_source == "absent"
        assert t.classify(0) is Callability.LOW_DEPTH
        assert t.classify(30) is Callability.HIGH_DEPTH

    def test_ceiling_excludes_a_doubled_locus(self):
        """The two-sided gate has to actually catch the thing it exists for: an
        unseparated paralogue doubles the depth at CN 2."""
        t = thresholds_for(2, LAMBDA1, DISPERSION)
        assert t.classify(2 * t.expected) is Callability.HIGH_DEPTH

    def test_alpha_widens_the_interval(self):
        tight = thresholds_for(2, LAMBDA1, DISPERSION, alpha=0.05)
        loose = thresholds_for(2, LAMBDA1, DISPERSION, alpha=0.0001)
        assert loose.ceiling > tight.ceiling
        assert loose.floor < tight.floor


class TestEffectiveCopies:
    def test_unique_sequence_is_the_genes_own_copies(self):
        assert effective_copies(2, shared_fraction=0.0, paralog_copies=3) == 2

    def test_fully_shared_block_counts_both_genes(self):
        assert effective_copies(2, shared_fraction=1.0, paralog_copies=3) == 5

    def test_partially_shared_interpolates(self):
        assert effective_copies(2, shared_fraction=0.5, paralog_copies=4) == 4

    def test_fraction_is_clamped(self):
        assert effective_copies(2, shared_fraction=1.7, paralog_copies=2) == 4
        assert effective_copies(2, shared_fraction=-0.3, paralog_copies=2) == 2


class TestCallPosition:
    def test_typical_diploid_position_is_callable(self):
        c = call_position(30, copies=2, lambda1=LAMBDA1, dispersion=DISPERSION)
        assert c.status is Callability.OK
        assert c.usable

    def test_low_depth_is_low_depth(self):
        c = call_position(3, copies=2, lambda1=LAMBDA1, dispersion=DISPERSION)
        assert c.status is Callability.LOW_DEPTH

    def test_pileup_at_unique_sequence_is_high_depth(self):
        c = call_position(120, copies=2, lambda1=LAMBDA1, dispersion=DISPERSION)
        assert c.status is Callability.HIGH_DEPTH

    def test_shared_block_is_not_rejected_for_being_deep(self):
        """LILRA6 CN 2 sharing a block with LILRB3 CN 2 sees ~4 copies' worth of
        reads. Judged against its own copy number that is a gross over-coverage;
        judged correctly it is exactly what should happen.

        Getting this wrong would reject ~90% of the LILRA6 coding sequence in
        every sample, for a modelling reason rather than a data one.
        """
        depth = 4 * LAMBDA1
        naive = call_position(depth, copies=2, lambda1=LAMBDA1, dispersion=DISPERSION)
        assert naive.status is Callability.HIGH_DEPTH

        aware = call_position(depth, copies=2, lambda1=LAMBDA1, dispersion=DISPERSION,
                              shared_fraction=1.0, paralog_copies=2)
        assert aware.status is Callability.OK
        assert aware.effective_copies == 4

    def test_shared_block_beyond_the_pair_is_flagged_as_ambiguous(self):
        """Over-covered even after accounting for the paralogue: assignment let
        in reads from somewhere the model does not know about."""
        c = call_position(300, copies=2, lambda1=LAMBDA1, dispersion=DISPERSION,
                          shared_fraction=1.0, paralog_copies=2)
        assert c.status is Callability.PARALOG_AMBIGUOUS

    def test_low_mapq_beats_depth(self):
        """A position whose reads nearly all multi-map is uninterpretable at any
        depth, and saying 'low depth' would send someone to look at coverage."""
        c = call_position(30, copies=2, lambda1=LAMBDA1, dispersion=DISPERSION,
                          mapq_fraction=0.1)
        assert c.status is Callability.LOW_MAPQ

    def test_no_coverage_model_is_not_a_zero(self):
        c = call_position(30, copies=2, lambda1=0.0)
        assert c.status is Callability.NO_MODEL
        assert not c.usable


class TestSummarise:
    def test_counts_and_fraction(self):
        calls = [
            call_position(30, 2, LAMBDA1, dispersion=DISPERSION),
            call_position(30, 2, LAMBDA1, dispersion=DISPERSION),
            call_position(2, 2, LAMBDA1, dispersion=DISPERSION),
            call_position(200, 2, LAMBDA1, dispersion=DISPERSION),
        ]
        s = summarise(calls)
        assert s["n_positions"] == 4
        assert s["n_ok"] == 2
        assert s["n_low_depth"] == 1
        assert s["n_high_depth"] == 1
        assert s["callable_fraction"] == 0.5

    def test_empty_is_zero_not_a_crash(self):
        assert summarise([])["callable_fraction"] == 0.0


class TestAgainstSimulation:
    """The model should keep ~(1 - 2*alpha) of positions that are behaving.

    Simulated from the same negative binomial the thresholds assume, so this
    checks the arithmetic and the normal approximation rather than the biology —
    whether real LILR depth is negative binomial is a question for validation
    against HPRC truth, not for a unit test.
    """

    @pytest.mark.parametrize("copies", [1, 2, 3, 4, 6])
    def test_retention_matches_alpha(self, copies):
        import random

        rng = random.Random(20260916)
        mean = copies * LAMBDA1
        r = DISPERSION * copies
        # Gamma-Poisson mixture: the standard construction of a negative binomial.
        depths = [rng.gammavariate(r, mean / r) for _ in range(20_000)]
        depths = [float(round(d)) for d in depths]

        calls = [call_position(d, copies, LAMBDA1, dispersion=DISPERSION) for d in depths]
        kept = summarise(calls)["callable_fraction"]

        expected = 1 - 2 * DEFAULT_ALPHA
        # The normal approximation is not exact in the tails at these means, and
        # the per-copy floor may bind; both move retention by a few percent.
        assert kept == pytest.approx(expected, abs=0.05), (
            f"copies={copies}: kept {kept:.3f}, expected ~{expected:.3f}"
        )
