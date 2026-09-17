"""Tests for the coverage model's pure functions.

The I/O path needs a BAM and is exercised in the integration tests; what is
checked here is the arithmetic that turns control-locus depth into λ₁, and the
diagnostic that decides whether MAPQ-20 numbers inside the LRC mean anything.

That diagnostic carries the most weight of anything in this file. If it is wrong
in the permissive direction, a non-ALT-aware cohort produces confident zeros at
LILRA6 and LILRB3 — which is a cohort of homozygous deletions, reported without
a warning, at genes where deletion is a real and interesting phenotype.
"""

from __future__ import annotations

import math
import statistics

import pytest

from lilrwgs.coverage import (
    MAX_DILUTION,
    MIN_BASES_PER_GC_BIN,
    MIN_Q20_LRC,
    ControlMeasurement,
    CoverageModel,
    alignment_verdict,
    fit_gc_correction,
    fit_lambda,
    gc_bin,
)


def control(name: str, mean_q20: float, *, mean_q0: float | None = None,
            var_q20: float | None = None, n_bases: int = 16_000,
            inside: bool = True) -> ControlMeasurement:
    """A control measurement with sensible defaults, Poisson-ish unless told."""
    return ControlMeasurement(
        name=name, chrom="chr19", start=0, end=n_bases,
        inside_placement=inside,
        mean_q0=mean_q0 if mean_q0 is not None else mean_q20,
        mean_q20=mean_q20,
        var_q20=var_q20 if var_q20 is not None else mean_q20,
        n_bases=n_bases,
    )


class TestFitLambda:
    def test_lambda_is_half_the_diploid_depth(self):
        controls = [control(f"C{i}", 30.0) for i in range(4)]
        lambda1, _ = fit_lambda(controls)
        assert lambda1 == pytest.approx(15.0)

    def test_tracks_library_depth(self):
        shallow, _ = fit_lambda([control("C", 20.0)])
        deep, _ = fit_lambda([control("C", 50.0)])
        assert shallow == pytest.approx(10.0)
        assert deep == pytest.approx(25.0)

    def test_median_resists_one_bad_control(self):
        """A control overlapping an unannotated CNV in one sample should move
        lambda by nothing; with four loci a mean would move it by a quarter of
        the error."""
        clean = [control(f"C{i}", 30.0) for i in range(3)]
        outlier = control("CNV", 60.0)
        lambda1, _ = fit_lambda(clean + [outlier])
        assert lambda1 == pytest.approx(15.0)

    def test_dispersion_is_halved_to_per_haploid(self):
        """NB is additive: a diploid locus is two haploid copies, so the
        per-haploid r is half the value fitted at the control's own depth."""
        # NB(mu=30, r=20) has variance 30 + 900/20 = 75.
        c = control("C", 30.0, var_q20=75.0, n_bases=100_000)
        _, dispersion = fit_lambda([c])
        assert dispersion == pytest.approx(10.0)

    def test_poisson_control_gives_infinite_dispersion(self):
        c = control("C", 30.0, var_q20=30.0, n_bases=100_000)
        _, dispersion = fit_lambda([c])
        assert dispersion == math.inf

    def test_no_usable_controls_is_zero_not_a_crash(self):
        lambda1, dispersion = fit_lambda([])
        assert lambda1 == 0.0
        assert dispersion == math.inf

    def test_zero_depth_controls_are_ignored(self):
        lambda1, _ = fit_lambda([control("dead", 0.0), control("live", 30.0)])
        assert lambda1 == pytest.approx(15.0)


class TestAlignmentVerdict:
    def test_alt_aware_is_the_good_case(self):
        assert alignment_verdict(q20_lrc=0.95, q20_outside=0.97,
                                 dilution=1.02) == "alt_aware"

    def test_lrc_collapsed_but_outside_alive_is_not_alt_aware(self):
        """The failure this exists to catch. MAPQ 20 dies inside the placement
        and survives outside it: nine alt haplotypes are taking the mapping
        quality, and every LILRA6/LILRB3 window reads as a deletion."""
        assert alignment_verdict(q20_lrc=0.02, q20_outside=0.96,
                                 dilution=1.1) == "not_alt_aware"

    def test_both_collapsed_is_a_different_problem(self):
        """If MAPQ is degraded outside the placement too, the alt contigs are
        not the cause and saying so would send someone to the wrong place."""
        assert alignment_verdict(q20_lrc=0.05, q20_outside=0.08,
                                 dilution=1.0) == "mapq_degraded_everywhere"

    def test_dilution_is_flagged_even_when_mapq_survives(self):
        assert alignment_verdict(q20_lrc=0.9, q20_outside=0.95,
                                 dilution=2.5) == "alt_diluted"

    @pytest.mark.parametrize("q20_lrc", [0.0, 0.1, MIN_Q20_LRC - 0.01])
    def test_threshold_is_respected(self, q20_lrc):
        assert alignment_verdict(q20_lrc, 0.95, 1.0) == "not_alt_aware"

    def test_just_above_threshold_passes(self):
        assert alignment_verdict(MIN_Q20_LRC + 0.01, 0.95, 1.0) == "alt_aware"


class TestUsableMapq20:
    def test_a_dead_lrc_is_reported_as_unusable_not_as_zero(self):
        """Dividing a dead MAPQ-20 count by a live baseline gives a clean,
        confident zero. A cohort of those looks exactly like a cohort of
        deletion homozygotes, so the model has to refuse rather than divide."""
        model = CoverageModel(sample="S", lambda1=15.0, q20_lrc=0.01,
                              q20_outside=0.96)
        assert not model.usable_mapq20

    def test_a_healthy_lrc_is_usable(self):
        model = CoverageModel(sample="S", lambda1=15.0, q20_lrc=0.94,
                              q20_outside=0.96)
        assert model.usable_mapq20


class TestGcCorrection:
    def test_flat_coverage_gives_flat_correction(self):
        n = MIN_BASES_PER_GC_BIN * 2
        gc = [0.40] * n + [0.50] * n
        depth = [30] * (2 * n)
        curve = fit_gc_correction(gc, depth)
        assert all(v == pytest.approx(1.0) for v in curve.values())

    def test_recovers_a_known_bias(self):
        """GC-rich positions covered at half depth should end up with half the
        multiplier of the GC-neutral ones, which is what stops λ₁ being
        systematically wrong at a GC-rich gene.

        The ratio between bins is the meaningful quantity; the absolute level
        depends on the normaliser, which is checked separately below.
        """
        n = MIN_BASES_PER_GC_BIN * 2
        gc = [0.40] * n + [0.65] * n
        depth = [30] * n + [15] * n
        curve = fit_gc_correction(gc, depth)
        assert curve[gc_bin(0.65)] / curve[gc_bin(0.40)] == pytest.approx(0.5, abs=0.02)

    def test_curve_composes_with_lambda(self):
        """λ₁ is the mean control depth over two, so applying the curve across
        the control GC distribution has to average back to λ₁.

        If the curve were normalised by the median instead, ``lambda_at`` would
        disagree with ``lambda1`` by the mean-median gap — a silent,
        GC-dependent bias in every threshold the pipeline sets.
        """
        n = MIN_BASES_PER_GC_BIN * 2
        gc = [0.40] * n + [0.65] * n
        depth = [30] * n + [15] * n
        curve = fit_gc_correction(gc, depth)

        model = CoverageModel(sample="S", lambda1=statistics.fmean(depth) / 2,
                              gc_correction=curve)
        recovered = statistics.fmean([model.lambda_at(g) for g in gc])
        assert recovered == pytest.approx(model.lambda1, rel=0.02)

    def test_thin_bins_are_omitted(self):
        """A noisy multiplier on λ₁ is worse than none: it would move thresholds
        for reasons unrelated to the sample."""
        n = MIN_BASES_PER_GC_BIN * 2
        gc = [0.40] * n + [0.80] * 10
        depth = [30] * n + [2] * 10
        curve = fit_gc_correction(gc, depth)
        assert gc_bin(0.80) not in curve

    def test_empty_input_is_empty_output(self):
        assert fit_gc_correction([], []) == {}


class TestLambdaAt:
    def test_no_curve_means_no_correction(self):
        model = CoverageModel(sample="S", lambda1=15.0)
        assert model.lambda_at(0.65) == 15.0

    def test_applies_the_bin_multiplier(self):
        model = CoverageModel(sample="S", lambda1=15.0,
                              gc_correction={gc_bin(0.65): 0.5})
        assert model.lambda_at(0.65) == pytest.approx(7.5)

    def test_unknown_bin_falls_back_to_uncorrected(self):
        model = CoverageModel(sample="S", lambda1=15.0,
                              gc_correction={gc_bin(0.65): 0.5})
        assert model.lambda_at(0.40) == 15.0


class TestGcBin:
    def test_bins_are_five_percent_wide(self):
        assert gc_bin(0.40) == gc_bin(0.41)
        assert gc_bin(0.40) != gc_bin(0.46)

    def test_clamped_to_range(self):
        assert 0 <= gc_bin(0.0) <= 20
        assert 0 <= gc_bin(1.0) <= 20
        assert 0 <= gc_bin(1.5) <= 20
