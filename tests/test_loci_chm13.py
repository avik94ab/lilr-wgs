"""Tests for the CHM13v2.0 coordinate basis.

Two properties carry the weight here.

The first is that CHM13 is a *different assembly with the same contig names*. A
GRCh38 coordinate used against CHM13 does not fail — it silently reads sequence
about 3 Mb away, which in this region is a different LILR gene. So the tests
assert that the two tables cannot be confused and that nothing in the CHM13 one
is a GRCh38 number left behind.

The second is that overlapping depth intervals double-count. That bug put λ₁ at
36.8 for a 30x library and halved every copy number in the cohort; it is the
reason `loci.merge_intervals` exists, and the derivation that produced these
windows hit it too — two runs of unique 100-mers separated by a gap shorter than
the window produce intervals that overlap by construction.
"""

from __future__ import annotations

import pytest

from lilrwgs import loci, loci_chm13

# Measured from chm13v2.0.fa.fai and the GRCh38 analysis set's .fai.
CHM13_CHR19 = 61_707_364
GRCH38_CHR19 = 58_617_616

# chr19:57,817,388 in CHM13, from the RefSeq Liftoff annotation. The slice stops
# short of it on purpose: KIR is the most copy-number-variable region in the
# genome and extracting it would multiply the readset for nothing.
KIR3DL3_START = 57_817_388


class TestNoOverlap:
    def test_slice_intervals_do_not_overlap(self):
        """The property whose violation doubled λ₁ and halved every copy number."""
        ivs = loci_chm13.slice_intervals()
        for a, b in zip(ivs, ivs[1:], strict=False):
            if a[0] == b[0]:
                assert a[2] <= b[1], f"{a} overlaps {b}"

    @pytest.mark.parametrize("gene", ["LILRA6", "LILRB3"])
    def test_unique_windows_do_not_overlap(self, gene):
        windows = sorted(loci_chm13.UNIQUE_WINDOWS[gene])
        for a, b in zip(windows, windows[1:], strict=False):
            assert a[2] <= b[1], f"{gene}: {a} overlaps {b}"

    @pytest.mark.parametrize("gene", ["LILRA6", "LILRB3"])
    def test_unique_windows_are_inside_their_gene(self, gene):
        _, g_start, g_end = loci_chm13.gene_span(gene)
        for _, s, e in loci_chm13.UNIQUE_WINDOWS[gene]:
            assert g_start <= s < e <= g_end


class TestNotGRCh38:
    """Nothing here may be a GRCh38 coordinate that was never updated."""

    def test_the_two_assemblies_are_named_apart(self):
        assert loci_chm13.ASSEMBLY != loci.ASSEMBLY

    def test_chr19_lengths_distinguish_them(self):
        assert loci_chm13.CHR19_LENGTH == CHM13_CHR19
        assert loci_chm13.CHR19_LENGTH != GRCH38_CHR19

    def test_no_gene_body_sits_at_its_grch38_address(self):
        for name, g in loci_chm13.GENE_BODIES.items():
            other = loci.GENE_BODIES[name]
            if other.on_primary:
                assert g.start != other.start, f"{name} kept its GRCh38 start"

    def test_the_slice_does_not_overlap_the_grch38_slice(self):
        """They are ~3 Mb apart; an overlap means one of them is wrong."""
        _, cs, ce = loci_chm13.LRC_SLICE
        _, gs, ge = loci.LRC_SLICE
        assert cs >= ge or ce <= gs

    def test_every_interval_fits_on_chm13_chr19(self):
        for _, _, end in loci_chm13.slice_intervals():
            assert end <= CHM13_CHR19


class TestLilra3IsOnThePrimaryAssembly:
    """The entire reason this module exists."""

    def test_lilra3_has_a_span(self):
        chrom, start, end = loci_chm13.gene_span("LILRA3")
        assert chrom == "chr19" and end > start

    def test_grch38_refuses_the_same_call(self):
        """The contrast that motivates the second assembly."""
        with pytest.raises(ValueError):
            loci.gene_span("LILRA3")

    def test_lilra3_length_matches_the_calling_reference(self):
        _, start, end = loci_chm13.gene_span("LILRA3")
        assert end - start == loci_chm13.LILRA3_SPAN == 7_126

    def test_lilra3_is_inside_the_extracted_slice(self):
        """If it were outside, every sample would read as a deletion homozygote
        — the same failure GRCh38's alt contigs cause when they are skipped."""
        _, s, e = loci_chm13.gene_span("LILRA3")
        _, ls, le = loci_chm13.LRC_SLICE
        assert ls <= s and e <= le

    def test_lilra3_sits_between_lilrb2_and_lilra5(self):
        """Its position in the cluster, which is how it was found."""
        _, _, b2_end = loci_chm13.gene_span("LILRB2")
        _, a5_start, _ = loci_chm13.gene_span("LILRA5")
        _, s, e = loci_chm13.gene_span("LILRA3")
        assert b2_end < s and e < a5_start


class TestSliceBounds:
    def test_kir_is_excluded(self):
        assert loci_chm13.LRC_SLICE[2] < KIR3DL3_START

    def test_every_lilr_gene_is_inside_the_slice(self):
        _, ls, le = loci_chm13.LRC_SLICE
        for name in loci_chm13.GENES:
            _, s, e = loci_chm13.gene_span(name)
            assert ls <= s and e <= le, f"{name} is outside the slice"

    def test_every_control_is_inside_some_slice_interval(self):
        """λ₁ comes from the controls; a slice without them has no baseline."""
        ivs = loci_chm13.slice_intervals()
        for c in loci_chm13.ALL_CONTROLS:
            assert any(s <= c.start and c.end <= e for _, s, e in ivs), \
                f"{c.name} is not covered"


class TestDropInForLoci:
    """`coverage` and `cn` take this module in place of `loci`, so the names
    they reach for have to exist and behave."""

    @pytest.mark.parametrize("name", [
        "CHROM", "GENES", "GENE_BODIES", "UNIQUE_WINDOWS", "CONTROL_LOCI",
        "OUTSIDE_LOCI", "ALL_CONTROLS", "MAPQ_STRICT", "MAPQ_ANY",
        "gene_span", "slice_intervals", "slice_regions", "extraction_regions",
        "as_region", "merge_intervals",
    ])
    def test_the_name_exists(self, name):
        assert hasattr(loci_chm13, name), f"coverage/cn reach for loci.{name}"

    def test_include_alts_is_accepted_and_ignored(self):
        """CHM13 has no alt contigs, but callers still pass the flag."""
        assert (loci_chm13.slice_intervals(include_alts=True)
                == loci_chm13.slice_intervals(include_alts=False))

    def test_mapq_floors_are_shared_with_grch38(self):
        """A threshold that differed between assemblies would make the two
        runs incomparable, which is the whole point of the coherence test."""
        assert loci_chm13.MAPQ_STRICT == loci.MAPQ_STRICT
        assert loci_chm13.MAPQ_ANY == loci.MAPQ_ANY


class TestUniqueWindowsMatchTheValidatedOnes:
    """The windows were transferred from GRCh38 by aligning their sequence, not
    re-derived, because re-deriving with the documented criterion does not
    reproduce the windows LILRA6's accuracy was validated on. Length is the
    check that the transfer moved the same sequence."""

    @pytest.mark.parametrize("gene", ["LILRA6", "LILRB3"])
    def test_total_length_is_preserved(self, gene):
        here = sum(e - s for _, s, e in loci_chm13.UNIQUE_WINDOWS[gene])
        there = sum(e - s for _, s, e in loci.UNIQUE_WINDOWS[gene])
        assert here == there

    @pytest.mark.parametrize("gene", ["LILRA6", "LILRB3"])
    def test_each_window_is_preserved(self, gene):
        here = [e - s for _, s, e in loci_chm13.UNIQUE_WINDOWS[gene]]
        there = [e - s for _, s, e in loci.UNIQUE_WINDOWS[gene]]
        assert here == there
