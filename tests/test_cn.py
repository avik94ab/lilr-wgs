"""Tests for copy-number calling.

The behaviour worth defending here is refusal. Two of these genes are measurable
only under conditions that do not always hold, and in both cases the failure mode
produces a confident wrong answer rather than an error: a non-ALT-aware alignment
turns LILRA6 into a homozygous deletion, and a failed region query turns LILRA3
into one. Since deletion is a real and common state at both genes, neither is
detectable downstream. So "not measured" has to survive as a distinct value.
"""

from __future__ import annotations

from lilrwgs.cn import (
    AMBIGUOUS_BAND,
    CNCall,
    CN_RANGE,
    _fit_unit,
    _call_unique_window,
    integerise,
    refine_cohort,
    tally_junction,
)
from lilrwgs.coverage import CoverageModel


def model(**kw) -> CoverageModel:
    base = dict(sample="S", lambda1=15.0, q20_lrc=0.95, q20_outside=0.96,
                dilution=1.0, alt_verdict="alt_aware")
    base.update(kw)
    return CoverageModel(**base)


class TestIntegerise:
    def test_rounds_to_the_nearest_copy(self):
        assert integerise(1.97, "LILRA6")[0] == 2
        assert integerise(3.04, "LILRA6")[0] == 3

    def test_confidence_is_highest_on_an_integer(self):
        assert integerise(2.0, "LILRA6")[1] == 1.0

    def test_confidence_is_lowest_on_a_boundary(self):
        assert integerise(2.5, "LILRA6")[1] == 0.0

    def test_clamped_to_the_observed_range(self):
        lo, hi = CN_RANGE["LILRA3"]
        assert integerise(9.0, "LILRA3")[0] == hi
        assert integerise(-2.0, "LILRA3")[0] == lo

    def test_out_of_range_is_not_confident(self):
        """Clamping must not turn an impossible estimate into a confident call
        at the boundary."""
        assert integerise(9.0, "LILRA3")[1] == 0.0

    def test_just_over_the_cap_is_still_confident(self):
        """2.047 at a gene capped at CN 2 rounds to 2 whether or not it is
        clamped, so it is a good call. Penalising it flagged four of five
        correct LILRA3 calls as ambiguous, and a flag that fires unevenly
        across copy-number classes biases any frequency computed from the
        confident subset."""
        copies, confidence = integerise(2.047, "LILRA3")
        assert copies == 2
        assert confidence > 0.85

    def test_ambiguous_band_is_flagged(self):
        call = CNCall(sample="S", gene="LILRA6", status="measured")
        call.estimate = 2.5
        call.copies, call.confidence = integerise(2.5, "LILRA6")
        assert call.ambiguous

        clean = CNCall(sample="S", gene="LILRA6", status="measured")
        clean.copies, clean.confidence = integerise(2.02, "LILRA6")
        assert not clean.ambiguous

    def test_unmeasured_is_never_ambiguous(self):
        """`ambiguous` qualifies a call. A gene that was not measured has no
        call to qualify, and reporting it as ambiguous would imply one exists."""
        call = CNCall(sample="S", gene="LILRA6", status="not_measured")
        assert not call.ambiguous


class TestRefusal:
    def test_dead_mapq_yields_not_measured_not_zero(self):
        """The single most important behaviour in this module.

        In a non-ALT-aware alignment every MAPQ-20 window in the LRC reads near
        zero. Divided by a live baseline that is a clean 0.0 copies, and a cohort
        of those is a cohort of LILRA6 deletion homozygotes — a wrong answer
        that arrives looking exactly like a right one.
        """
        call = _call_unique_window(
            "S", "LILRA6", bam="/nonexistent.bam",
            model=model(q20_lrc=0.01, alt_verdict="not_alt_aware"),
            reference=None, samtools="samtools",
        )
        assert call.status == "not_measured"
        assert call.copies is None
        assert call.estimate is None
        assert any("not measured rather than zero" in n for n in call.notes)

    def test_no_coverage_model_is_failed_not_zero(self):
        call = _call_unique_window(
            "S", "LILRA6", bam="/nonexistent.bam", model=model(lambda1=0.0),
            reference=None, samtools="samtools",
        )
        assert call.status == "failed"
        assert call.copies is None

    def test_failed_query_is_failed_not_zero(self):
        """A region query that fails and a region that is genuinely empty are
        different facts. At LILRA3 the second is a common true answer."""
        call = _call_unique_window(
            "S", "LILRA6", bam="/definitely/not/a.bam", model=model(),
            reference=None, samtools="samtools",
        )
        assert call.status == "failed"
        assert call.copies is None


class TestCohortRefinement:
    def _calls(self, gene: str, estimates: list[float]) -> list[CNCall]:
        out = []
        for i, e in enumerate(estimates):
            c = CNCall(sample=f"S{i}", gene=gene, status="measured", estimate=e)
            c.copies, c.confidence = integerise(e, gene)
            out.append(c)
        return out

    def test_reports_a_correct_scale_as_correct(self):
        estimates = [2.0, 2.05, 1.95, 3.0, 3.02, 4.0, 2.0, 2.1] * 4
        result = refine_cohort(self._calls("LILRA6", estimates))
        assert result["LILRA6"]["unit"] is not None
        assert abs(result["LILRA6"]["unit"] - 1.0) < 0.08
        assert "scale looks right" in result["LILRA6"]["note"]

    def test_detects_a_systematic_offset(self):
        """If every estimate is 20% high, lambda1 is wrong by that factor and
        the cohort should say so rather than quietly rescale."""
        estimates = [e * 1.2 for e in [2.0, 2.0, 3.0, 3.0, 4.0, 2.0, 2.0, 3.0] * 4]
        result = refine_cohort(self._calls("LILRA6", estimates))
        assert abs(result["LILRA6"]["unit"] - 1.2) < 0.1
        assert "systematically off" in result["LILRA6"]["note"]

    def test_small_cohort_declines_to_fit(self):
        result = refine_cohort(self._calls("LILRA6", [2.0, 2.0, 3.0]))
        assert result["LILRA6"]["unit"] is None
        assert "too few" in result["LILRA6"]["note"]

    def test_unmeasured_calls_are_excluded(self):
        calls = self._calls("LILRA6", [2.0] * 25)
        calls += [CNCall(sample="bad", gene="LILRA6", status="not_measured")]
        result = refine_cohort(calls)
        assert result["LILRA6"]["n"] == 25

    def test_a_single_sample_needs_no_cohort(self):
        """The structural claim of the rewrite: one sample is callable on its
        own, so refinement being unavailable must not be an error."""
        result = refine_cohort(self._calls("LILRA6", [2.03]))
        assert result["LILRA6"]["unit"] is None


class TestFitUnit:
    def test_recovers_unit_one(self):
        assert abs(_fit_unit([1.0, 2.0, 3.0, 2.0, 4.0, 2.0]) - 1.0) < 0.02

    def test_recovers_a_stretched_unit(self):
        assert abs(_fit_unit([1.15, 2.30, 3.45, 2.30, 4.60]) - 1.15) < 0.03


class TestDepthCommand:
    """The samtools command must be well-formed.

    A flag spliced in at the wrong index landed before the `depth` subcommand.
    samtools rejected it, `_mean_depth` turned the failure into None, and the
    symptom was every LILRA6 and LILRB3 call coming back "failed" — nothing in
    the output pointed at a malformed command line.
    """

    def _command(self, monkeypatch, **kwargs) -> list[str]:
        from lilrwgs import cn as cn_mod

        seen: list[list[str]] = []

        def fake_run(cmd, **_):
            seen.append(cmd)
            class R:
                stdout = "chr19\t100\t30\n"
            return R()

        monkeypatch.setattr(cn_mod, "run", fake_run)
        monkeypatch.setattr(cn_mod, "require", lambda *a, **k: None)
        cn_mod._mean_depth("x.bam", [("chr19", 0, 10)], 20, **kwargs)
        return seen[0]

    def test_subcommand_comes_first(self, monkeypatch):
        cmd = self._command(monkeypatch)
        assert cmd[1] == "depth", f"argv[1] must be the subcommand, got {cmd[:4]}"

    def test_supplementary_excluded_by_default(self, monkeypatch):
        cmd = self._command(monkeypatch)
        assert "-G" in cmd and cmd[cmd.index("-G") + 1] == "0x800"

    def test_supplementary_kept_when_asked(self, monkeypatch):
        """LILRA3's evidence is almost entirely supplementary records, because
        bwa -Y puts an ALT-contig hit there and LILRA3 has no primary locus."""
        cmd = self._command(monkeypatch, supplementary=True)
        assert "-G" not in cmd
        assert cmd[1] == "depth"


class TestJunctionTally:
    """The LILRA3 deletion breakpoint.

    This assay silently measured nothing for the whole of the pilot. The recorded
    junction base was 28 bp left of where reads actually clip, so `clipped` came
    back 0–4 whatever the donor's copy number was, `LILRA3` was called from depth
    alone, and — because the two routes then disagreed by construction — every
    LILRA3-bearing donor picked up a "treat this call as unresolved" note. The
    numbers below are from the 101-donor HPRC overlap at the corrected base.
    """

    LEFT = 54_297_005
    RIGHT = LEFT + 5          # the microhomology's far end

    def sam(self, pos: int, cigar: str) -> str:
        """A minimal SAM record; only POS and CIGAR are read."""
        return f"r\t0\tchr19\t{pos}\t60\t{cigar}\t*\t0\t0\t*\t*"

    def test_read_ending_clipped_at_the_breakpoint_counts_as_clipped(self):
        # 100 aligned bases arriving at LEFT, then clipped into LILRA3.
        lines = [self.sam(self.LEFT - 100 + 1, "100M50S")]
        assert tally_junction(lines, self.LEFT, self.RIGHT) == (1, 0)

    def test_read_starting_clipped_at_the_breakpoint_counts_as_clipped(self):
        # Coming back out of LILRA3 onto the right flank.
        lines = [self.sam(self.RIGHT + 1, "50S100M")]
        assert tally_junction(lines, self.LEFT, self.RIGHT) == (1, 0)

    def test_read_crossing_cleanly_counts_as_spanning(self):
        lines = [self.sam(self.LEFT - 75 + 1, "150M")]
        assert tally_junction(lines, self.LEFT, self.RIGHT) == (0, 1)

    def test_the_old_junction_base_is_not_where_reads_clip(self):
        """The regression itself: a clip at the old base is not the signal.

        Reads clipping 28 bp away are ordinary alignment noise, and counting them
        is what made the assay look alive while measuring nothing.
        """
        old = 54_296_977
        lines = [self.sam(old - 100 + 1, "100M50S")]
        assert tally_junction(lines, self.LEFT, self.RIGHT) == (0, 0)

    def test_short_clips_are_not_the_signal(self):
        """Two or three trimmed bases are adapter, not a breakpoint."""
        lines = [self.sam(self.LEFT - 100 + 1, "100M4S")]
        clipped, _ = tally_junction(lines, self.LEFT, self.RIGHT)
        assert clipped == 0

    def test_a_clipped_read_is_not_also_counted_as_spanning(self):
        """The classes have to partition, or the ratio stops being a fraction."""
        lines = [self.sam(self.LEFT - 100 + 1, "100M50S")] * 3
        clipped, spanning = tally_junction(lines, self.LEFT, self.RIGHT)
        assert (clipped, spanning) == (3, 0)

    def test_deleted_homozygote_shape(self):
        """Truth CN 0: reads cross, none clip. Estimate 0.0."""
        lines = [self.sam(self.LEFT - 75 + 1, "150M")] * 30
        clipped, spanning = tally_junction(lines, self.LEFT, self.RIGHT)
        assert clipped == 0
        assert 2.0 * clipped / (clipped + spanning) == 0.0

    def test_bearing_homozygote_shape(self):
        """Truth CN 2: nothing crosses. Measured in 65/65 donors of the overlap."""
        lines = ([self.sam(self.LEFT - 100 + 1, "100M50S")] * 38
                 + [self.sam(self.RIGHT + 1, "50S100M")] * 19)
        clipped, spanning = tally_junction(lines, self.LEFT, self.RIGHT)
        assert spanning == 0
        assert 2.0 * clipped / (clipped + spanning) == 2.0

    def test_heterozygote_rounds_to_one(self):
        """Truth CN 1 reads a median 1.12 over the overlap — inside the band.

        The 12% excess is the bearing chromosome offering two breakpoints' worth
        of clipped reads against the deleted one's single spanning window.
        """
        lines = ([self.sam(self.LEFT - 100 + 1, "100M50S")] * 14
                 + [self.sam(self.RIGHT + 1, "50S100M")] * 8
                 + [self.sam(self.LEFT - 75 + 1, "150M")] * 18)
        clipped, spanning = tally_junction(lines, self.LEFT, self.RIGHT)
        assert (clipped, spanning) == (22, 18)     # HG00253, measured
        assert integerise(2.0 * clipped / (clipped + spanning), "LILRA3")[0] == 1


class TestUnitFitDeclinesWhenUnconstrained:
    """The unit is a spacing, so one copy-number class does not constrain it.

    LILRB3 over the 101-donor overlap is 99 donors at CN 2 and 2 at CN 1. The
    grid fitted 0.700 and reported lambda1 as "systematically off by that
    factor" — on the gene that scored 100% against truth.
    """

    def _calls(self, gene: str, estimates: list[float]) -> list[CNCall]:
        out = []
        for i, e in enumerate(estimates):
            c = CNCall(sample=f"S{i}", gene=gene, status="measured", estimate=e)
            c.copies, c.confidence = integerise(e, gene)
            out.append(c)
        return out

    def test_a_uniform_cohort_gets_no_unit(self):
        estimates = [2.06] * 99 + [1.03, 1.05]
        result = refine_cohort(self._calls("LILRB3", estimates))
        assert result["LILRB3"]["unit"] is None
        assert "too uniform" in result["LILRB3"]["note"]

    def test_a_uniform_cohort_still_reports_its_median(self):
        """Declining to fit is not declining to report. The median is what says
        the scale is fine, and it is the number a human would look at next."""
        result = refine_cohort(self._calls("LILRB3", [2.06] * 99 + [1.03, 1.05]))
        assert result["LILRB3"]["median_estimate"] == 2.06

    def test_two_populated_classes_are_enough_to_fit(self):
        """LILRA3 in the same cohort: 21 donors at CN 1, 67 at CN 2. It fits."""
        estimates = [1.02] * 21 + [1.98] * 67
        result = refine_cohort(self._calls("LILRA3", estimates))
        assert result["LILRA3"]["unit"] is not None
        assert abs(result["LILRA3"]["unit"] - 1.0) < 0.08

    def test_a_real_offset_is_still_caught(self):
        """The check must keep working where it is meaningful."""
        estimates = [e * 1.2 for e in [1.0] * 20 + [2.0] * 40 + [3.0] * 20]
        result = refine_cohort(self._calls("LILRA6", estimates))
        assert result["LILRA6"]["unit"] is not None
        assert "systematically off" in result["LILRA6"]["note"]
