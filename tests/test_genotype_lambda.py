"""λ₁ reaches the callable track in the right units, and the two copies agree.

`docs/variant_calling.md` §2.3. λ₁ is measured on the CRAM slice and applied to a
realigned per-gene BAM, after extraction, panel recruitment, arbitration and
realignment have each dropped reads. `recruitment_efficiency()` converts between
the two — 0.813 on HG00096, so roughly a fifth.

`CoverageModel.lambda_at()` owns that conversion. `genotype.main()` cannot call
it, because it reads the coverage JSON rather than holding the object, so it
reconstructs the same arithmetic inline. That is reasonable and it is also a
second copy, and the failure mode if the copies drift is invisible: every
threshold in the callable track moves by ~19% and no field in any output
changes. A gene simply reads as under-covered by the pipeline's own losses,
which is the predecessor's failure mode arriving by a different route.

These tests pin the conversion and assert the two copies agree.
"""

from __future__ import annotations

import inspect

import pytest

from lilrwgs import genotype
from lilrwgs.coverage import CoverageModel

# Measured on HG00096 and quoted throughout the docs.
HG00096_EFFICIENCY = 0.813
HG00096_LAMBDA1 = 19.1


def _model(lambda1=HG00096_LAMBDA1, efficiency=HG00096_EFFICIENCY,
           gc_correction=None):
    m = CoverageModel(sample="T")
    m.lambda1 = lambda1
    m.efficiency = efficiency
    m.gc_correction = gc_correction or {}
    return m


def _genotype_lambda(cov: dict) -> float:
    """The conversion exactly as `genotype.main()` performs it.

    Kept as a transcription rather than a call, because the point is to detect
    the inline copy changing. If this stops matching the source, the assertion
    in `test_the_inline_copy_is_still_the_same_expression` fails and says so.
    """
    efficiency = float(cov.get("efficiency", 1.0))
    return float(cov["lambda1"]) * efficiency


class TestTheConversionIsApplied:
    def test_efficiency_multiplies_lambda1(self):
        cov = {"lambda1": HG00096_LAMBDA1, "efficiency": HG00096_EFFICIENCY}
        assert _genotype_lambda(cov) == pytest.approx(15.53, abs=0.01)

    def test_dropping_it_would_move_thresholds_by_about_a_fifth(self):
        """The size of the error being guarded against, stated as a number."""
        cov = {"lambda1": HG00096_LAMBDA1, "efficiency": HG00096_EFFICIENCY}
        with_it = _genotype_lambda(cov)
        without = float(cov["lambda1"])
        assert (without - with_it) / without == pytest.approx(0.187, abs=0.005)

    def test_a_missing_efficiency_defaults_to_one_rather_than_zero(self):
        """An older coverage JSON has no `efficiency` key. Defaulting to 0 would
        make lambda_1 zero and every position NO_MODEL; defaulting to 1 leaves
        the thresholds merely strict, which is the safe direction."""
        assert _genotype_lambda({"lambda1": 18.0}) == pytest.approx(18.0)


class TestTheTwoCopiesAgree:
    """`CoverageModel.lambda_at()` and `genotype.main()`'s inline arithmetic
    must produce the same number, or the callable track is judged against a
    different λ₁ than everything else in the pipeline."""

    def test_they_agree_without_a_gc_correction(self):
        model = _model()
        cov = {"lambda1": model.lambda1, "efficiency": model.efficiency}
        assert _genotype_lambda(cov) == pytest.approx(model.lambda_at())

    @pytest.mark.parametrize("gc,factor", [(0.35, 0.9), (0.50, 1.0), (0.65, 1.1)])
    def test_they_agree_with_a_gc_correction(self, gc, factor):
        from lilrwgs.coverage import gc_bin

        correction = {gc_bin(gc): factor}
        model = _model(gc_correction=correction)
        cov = {"lambda1": model.lambda1, "efficiency": model.efficiency}

        # genotype.main()'s closure: the converted lambda_1, scaled by GC.
        inline = _genotype_lambda(cov) * correction.get(gc_bin(gc), 1.0)
        assert inline == pytest.approx(model.lambda_at(gc))

    def test_an_unseen_gc_bin_falls_back_to_one_in_both(self):
        """A bin with too little support is omitted from the correction, and a
        position landing in it must be judged at the uncorrected λ₁ rather than
        at zero."""
        from lilrwgs.coverage import gc_bin

        model = _model(gc_correction={gc_bin(0.5): 1.2})
        cov = {"lambda1": model.lambda1, "efficiency": model.efficiency}
        inline = _genotype_lambda(cov) * model.gc_correction.get(gc_bin(0.9), 1.0)
        assert inline == pytest.approx(model.lambda_at(0.9))
        assert inline == pytest.approx(model.lambda1 * model.efficiency)


class TestTheInlineCopyHasNotDrifted:
    """A transcription test. If `genotype.main()` stops computing λ₁ this way,
    these tests are asserting something the code no longer does, and the failure
    should point at the transcription rather than at the arithmetic."""

    def test_the_inline_copy_is_still_the_same_expression(self):
        src = inspect.getsource(genotype.main)
        assert 'efficiency = float(cov.get("efficiency", 1.0))' in src
        assert 'lambda1 = float(cov["lambda1"]) * efficiency' in src

    def test_the_gc_closure_scales_the_converted_value(self):
        """Not the raw lambda1 -- scaling the unconverted value would apply the
        GC correction and drop the efficiency factor in one step."""
        src = inspect.getsource(genotype.main)
        assert "return lambda1 * gc.get(gc_bin(value), 1.0)" in src

    def test_the_reason_is_recorded_next_to_the_code(self):
        """CLAUDE.md's convention: the measurement lives beside the constant it
        justifies, not only in a commit message."""
        src = inspect.getsource(genotype.main)
        assert "realigned per-gene BAM" in src
