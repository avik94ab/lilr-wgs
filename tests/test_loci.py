"""Tests for the coordinate basis.

The interval-merging tests exist because of a bug that produced entirely
plausible output. `samtools view` given several regions emits a read once per
region it overlaps, and three control loci sit inside the LRC slice — so their
depth doubled, lambda_1 came out at 36.8 for a 30x library, and every copy
number derived from it was halved. A cohort of LILRA6 hemizygotes is a
believable result, which is exactly why nothing downstream would have caught it.
"""

from __future__ import annotations

from lilrwgs import loci


class TestMergeIntervals:
    def test_contained_interval_is_absorbed(self):
        merged = loci.merge_intervals([("chr19", 100, 1000), ("chr19", 200, 300)])
        assert merged == [("chr19", 100, 1000)]

    def test_overlapping_intervals_are_joined(self):
        merged = loci.merge_intervals([("chr19", 100, 300), ("chr19", 250, 500)])
        assert merged == [("chr19", 100, 500)]

    def test_touching_intervals_are_joined(self):
        merged = loci.merge_intervals([("chr19", 100, 300), ("chr19", 300, 500)])
        assert merged == [("chr19", 100, 500)]

    def test_disjoint_intervals_are_kept_apart(self):
        merged = loci.merge_intervals([("chr19", 100, 200), ("chr19", 400, 500)])
        assert merged == [("chr19", 100, 200), ("chr19", 400, 500)]

    def test_contigs_do_not_merge_across(self):
        merged = loci.merge_intervals([("chr19", 100, 200), ("chr1", 150, 250)])
        assert len(merged) == 2

    def test_unsorted_input(self):
        merged = loci.merge_intervals([("chr19", 400, 500), ("chr19", 100, 450)])
        assert merged == [("chr19", 100, 500)]


class TestSliceIntervals:
    def test_no_two_intervals_overlap(self):
        """The property that was violated in the real bug."""
        intervals = loci.slice_intervals()
        by_chrom: dict[str, list] = {}
        for chrom, start, end in intervals:
            by_chrom.setdefault(chrom, []).append((start, end))
        for chrom, spans in by_chrom.items():
            spans.sort()
            for (_, prev_end), (next_start, _) in zip(spans, spans[1:]):
                assert next_start >= prev_end, f"{chrom} intervals overlap"

    def test_every_control_is_covered(self):
        """lambda_1 is meaningless if a control is not in the slice."""
        intervals = loci.slice_intervals()
        for control in loci.ALL_CONTROLS:
            assert any(c == control.chrom and s <= control.start and e >= control.end
                       for c, s, e in intervals), f"{control.name} is not in the slice"

    def test_the_junction_is_covered(self):
        chrom, pos = loci.LILRA3_JUNCTION
        assert any(c == chrom and s <= pos < e for c, s, e in loci.slice_intervals())

    def test_alt_contigs_are_included_by_default(self):
        contigs = {c for c, _, _ in loci.slice_intervals()}
        assert any(c.endswith("_alt") for c in contigs)

    def test_alt_contigs_can_be_excluded(self):
        contigs = {c for c, _, _ in loci.slice_intervals(include_alts=False)}
        assert not any(c.endswith("_alt") for c in contigs)


class TestGeneBodies:
    def test_lilra3_has_no_primary_span(self):
        """It is deleted from the GRCh38 primary assembly, and code that reaches
        gene_span() with it has made an assumption that does not hold."""
        import pytest
        with pytest.raises(ValueError, match="primary assembly"):
            loci.gene_span("LILRA3")

    def test_unique_windows_sit_inside_their_gene(self):
        for gene, windows in loci.UNIQUE_WINDOWS.items():
            chrom, start, end = loci.gene_span(gene)
            for c, s, e in windows:
                assert c == chrom and s >= start and e <= end

    def test_all_eleven_genes_are_present(self):
        assert len(loci.GENES) == 11
        assert set(loci.GENES) == set(loci.GENE_BODIES)
