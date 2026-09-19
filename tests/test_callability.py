"""Tests for the per-position evidence track.

This module is the junction between arithmetic that is verified and I/O that was
not. `depth_model` is pure and well covered; `callability` is what feeds it real
numbers off a BAM, and every one of its decisions could be individually wrong
without the output looking wrong — a depth that silently disagrees with
`samtools depth`, a base-quality filter biting harder than intended, a
shared-read count that misses reads whose tag did not survive the FASTQ hop.

The most valuable test here is `TestDepthAgreesWithSamtools`, because a wrong
depth produces a wrong callable track and a plausible consensus, and nothing
downstream can tell.

`docs/variant_calling.md` §1.1, §2.1, §2.2.
"""

from __future__ import annotations

import subprocess

import pytest

from lilrwgs import callability
from lilrwgs.assign import SHARED_TAG
from lilrwgs.callability import PositionEvidence
from lilrwgs.depth_model import Callability, call_position

pysam = pytest.importorskip("pysam")

REF_LEN = 300
CONTIG = "LILRB1"


# --------------------------------------------------------------------------
# fixtures: small BAMs with known, hand-countable content
# --------------------------------------------------------------------------

def _header():
    return {"HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"LN": REF_LEN, "SN": CONTIG}]}


def _read(name, pos, length=50, mapq=60, flags=0, cigar=None, tags=None,
          baseq=30):
    a = pysam.AlignedSegment()
    a.query_name = name
    a.query_sequence = "A" * length
    a.flag = flags
    a.reference_id = 0
    a.reference_start = pos
    a.mapping_quality = mapq
    a.cigar = cigar or [(0, length)]          # 0 = M
    a.query_qualities = pysam.qualitystring_to_array(chr(33 + baseq) * length)
    for tag, value in (tags or {}).items():
        a.set_tag(tag, value)
    return a


def _write_bam(path, reads):
    with pysam.AlignmentFile(str(path), "wb", header=_header()) as out:
        for r in reads:
            out.write(r)
    pysam.index(str(path))
    return str(path)


@pytest.fixture
def simple_bam(tmp_path):
    """Ten 50 bp reads stacked at position 100, all clean."""
    reads = [_read(f"r{i}", 100) for i in range(10)]
    return _write_bam(tmp_path / "simple.bam", reads)


# --------------------------------------------------------------------------


class TestDepthAgreesWithSamtools:
    """The junction test. `gather_evidence` must count what `samtools depth`
    counts, or the whole model is being applied to the wrong numbers."""

    def _samtools_depth(self, bam, extra=()):
        out = subprocess.run(
            ["samtools", "depth", "-a", "-q", "13", "-G", "0x900",
             *extra, bam],
            capture_output=True, text=True, check=True).stdout
        depths = {}
        for line in out.splitlines():
            f = line.split("\t")
            if len(f) >= 3:
                depths[int(f[1]) - 1] = int(f[2])    # to 0-based
        return depths

    def test_flat_pileup_matches_position_by_position(self, simple_bam):
        ours = callability.gather_evidence(simple_bam, REF_LEN)
        theirs = self._samtools_depth(simple_bam)
        mismatches = [(e.pos, e.depth, theirs.get(e.pos, 0))
                      for e in ours if e.depth != theirs.get(e.pos, 0)]
        assert not mismatches, f"first 5 disagreements: {mismatches[:5]}"

    def test_ragged_pileup_matches(self, tmp_path):
        """Staggered starts, so every edge is exercised rather than one block."""
        reads = [_read(f"r{i}", 50 + 7 * i, length=40) for i in range(12)]
        bam = _write_bam(tmp_path / "ragged.bam", reads)
        ours = callability.gather_evidence(bam, REF_LEN)
        theirs = self._samtools_depth(bam)
        bad = [(e.pos, e.depth, theirs.get(e.pos, 0))
               for e in ours if e.depth != theirs.get(e.pos, 0)]
        assert not bad, f"first 5: {bad[:5]}"

    def test_deletions_are_not_counted_as_depth(self, tmp_path):
        """A read spanning a deletion covers the flanks and not the gap.

        `is_del` is skipped in gather_evidence; samtools does the same, and the
        two must agree or a deleted stretch reads as covered.
        """
        cig = [(0, 20), (2, 10), (0, 20)]          # 20M 10D 20M
        reads = [_read(f"d{i}", 100, length=40, cigar=cig) for i in range(8)]
        bam = _write_bam(tmp_path / "del.bam", reads)
        ours = callability.gather_evidence(bam, REF_LEN)
        theirs = self._samtools_depth(bam)
        for p in range(120, 130):                   # inside the deletion
            assert ours[p].depth == 0, f"pos {p} counted a deleted base"
        bad = [(e.pos, e.depth, theirs.get(e.pos, 0))
               for e in ours if e.depth != theirs.get(e.pos, 0)]
        assert not bad, f"first 5: {bad[:5]}"

    def test_refskips_are_not_counted(self, tmp_path):
        cig = [(0, 20), (3, 30), (0, 20)]          # 20M 30N 20M
        reads = [_read(f"n{i}", 100, length=40, cigar=cig) for i in range(6)]
        bam = _write_bam(tmp_path / "skip.bam", reads)
        ours = callability.gather_evidence(bam, REF_LEN)
        for p in range(120, 150):
            assert ours[p].depth == 0, f"pos {p} counted a refskip"

    def test_secondary_and_supplementary_are_excluded(self, tmp_path):
        """Both are fragments of a read already counted; counting them again
        inflates depth exactly where paralogues make alignments ambiguous."""
        reads = [_read("primary", 100)]
        reads.append(_read("sec", 100, flags=256))          # secondary
        reads.append(_read("sup", 100, flags=2048))         # supplementary
        bam = _write_bam(tmp_path / "flags.bam", reads)
        ours = callability.gather_evidence(bam, REF_LEN)
        assert ours[110].depth == 1


class TestMapqAndSharedCounting:
    def test_mapq_floor_is_counted_separately_from_depth(self, tmp_path):
        """A low-MAPQ read still contributes depth; it just does not vouch for
        it. Conflating the two hides multi-mapping as missing coverage."""
        reads = [_read("hi", 100, mapq=60), _read("lo", 100, mapq=5)]
        bam = _write_bam(tmp_path / "mapq.bam", reads)
        e = callability.gather_evidence(bam, REF_LEN)[110]
        assert e.depth == 2
        assert e.n_pass_mapq == 1
        assert e.mapq_fraction == pytest.approx(0.5)

    def test_shared_reads_are_found_by_tag(self, tmp_path):
        reads = [_read("plain", 100),
                 _read("shared", 100, tags={SHARED_TAG: 2})]
        bam = _write_bam(tmp_path / "zs.bam", reads)
        e = callability.gather_evidence(bam, REF_LEN)[110]
        assert e.n_shared == 1

    def test_shared_reads_are_found_by_name_after_the_fastq_hop(self, tmp_path):
        """The ZS tag does not survive BAM -> FASTQ -> realignment, which is
        exactly the path this pipeline takes. `shared_names` is how the fact
        crosses, and if it did not work the shared fraction would read 0 and the
        two-sided ceiling would reject real shared-block positions as pile-ups.
        """
        reads = [_read("plain", 100), _read("was_shared", 100)]
        bam = _write_bam(tmp_path / "names.bam", reads)
        e = callability.gather_evidence(
            bam, REF_LEN, shared_names={"was_shared": 2})[110]
        assert e.n_shared == 1

    def test_tag_and_table_do_not_double_count(self, tmp_path):
        """A read both tagged and listed is one read, not two."""
        reads = [_read("both", 100, tags={SHARED_TAG: 2})]
        bam = _write_bam(tmp_path / "both.bam", reads)
        e = callability.gather_evidence(
            bam, REF_LEN, shared_names={"both": 2})[110]
        assert e.depth == 1 and e.n_shared == 1

    def test_base_quality_filter_applies(self, tmp_path):
        """min_baseq=13 matches what coverage.py measures lambda_1 with, so the
        threshold and the number it is compared against see the same bases."""
        reads = [_read("good", 100, baseq=30), _read("bad", 100, baseq=2)]
        bam = _write_bam(tmp_path / "baseq.bam", reads)
        e = callability.gather_evidence(bam, REF_LEN)[110]
        assert e.depth == 1


class TestEvidenceFractions:
    """Zero depth must not divide by zero; it must report 0.0 and let the model
    call it low_depth rather than raising."""

    def test_no_reads_gives_zero_fractions(self):
        e = PositionEvidence(pos=0)
        assert e.mapq_fraction == 0.0
        assert e.shared_fraction == 0.0

    def test_fractions_are_over_depth(self):
        e = PositionEvidence(pos=0, depth=8, n_pass_mapq=6, n_shared=2)
        assert e.mapq_fraction == pytest.approx(0.75)
        assert e.shared_fraction == pytest.approx(0.25)


class TestClassify:
    """`classify` is a thin router into depth_model; what it owns is applying
    the GC lookup per position rather than per gene."""

    def _evidence(self, depths):
        return [PositionEvidence(pos=i, depth=d, n_pass_mapq=d)
                for i, d in enumerate(depths)]

    def test_depth_inside_the_band_is_ok(self):
        calls = callability.classify(self._evidence([36] * 5), copies=2,
                                     lambda1=18.0, dispersion=float("inf"))
        assert all(c.status is Callability.OK for c in calls)

    def test_no_reads_at_all_is_low_depth_not_low_mapq(self):
        """At zero depth `mapq_fraction` is 0/0, reported as 0.0. Without a
        guard that fails the MAPQ test and an uncovered position is labelled
        LOW_MAPQ -- so a gene at copy number 0, a LILRA3 deletion homozygote,
        reports its whole length as a mapping problem. The two statuses
        prescribe different actions, which is the only reason to have both.
        """
        calls = callability.classify(self._evidence([0]), copies=2,
                                     lambda1=18.0, dispersion=float("inf"),
                                     min_mapq_fraction=0.5)
        assert calls[0].status is Callability.LOW_DEPTH

    def test_too_little_depth_is_low_depth_not_ok(self):
        calls = callability.classify(self._evidence([2]), copies=2,
                                     lambda1=18.0, dispersion=float("inf"))
        assert calls[0].status is Callability.LOW_DEPTH

    def test_too_much_depth_is_high_depth(self):
        """The two-sided gate: a pile-up is not better evidence."""
        calls = callability.classify(self._evidence([200]), copies=2,
                                     lambda1=18.0, dispersion=float("inf"))
        assert calls[0].status is Callability.HIGH_DEPTH

    def test_low_mapq_is_distinguished_from_low_depth(self):
        ev = [PositionEvidence(pos=0, depth=36, n_pass_mapq=2)]
        calls = callability.classify(ev, copies=2, lambda1=18.0,
                                     dispersion=float("inf"))
        assert calls[0].status is Callability.LOW_MAPQ

    def test_gc_lookup_is_applied_per_position(self):
        """A GC-extreme exon must be judged against what its own GC predicts,
        not against the gene's average — the whole reason the lookup exists."""
        seen = []

        def lookup(gc):
            seen.append(gc)
            return 18.0

        callability.classify(self._evidence([36, 36, 36]), copies=2,
                             lambda1=18.0, dispersion=float("inf"),
                             gc_by_pos=[0.3, 0.5, 0.7], gc_lookup=lookup)
        assert seen == [0.3, 0.5, 0.7]

    def test_without_a_lookup_the_flat_lambda_is_used(self):
        calls = callability.classify(self._evidence([36]), copies=2,
                                     lambda1=18.0, dispersion=float("inf"))
        expected = call_position(36, copies=2, lambda1=18.0,
                                 dispersion=float("inf"))
        assert calls[0].expected == expected.expected


class TestCallableBed:
    """The BED is what constrains HaplotypeCaller, so its arithmetic has to be
    right in the same way a coordinate does."""

    def _calls(self, statuses):
        return [call_position(36 if s else 0, copies=2, lambda1=18.0,
                              dispersion=float("inf")) for s in statuses]

    def test_runs_are_merged_and_half_open(self, tmp_path):
        calls = self._calls([1, 1, 1, 0, 0, 1, 1])
        bed = tmp_path / "c.bed"
        n = callability.callable_bed(bed, CONTIG, calls)
        rows = [ln.split("\t") for ln in bed.read_text().splitlines()]
        assert [(int(r[1]), int(r[2])) for r in rows] == [(0, 3), (5, 7)]
        assert n == 5

    def test_a_trailing_run_is_closed(self, tmp_path):
        """An off-by-one here silently drops the last callable stretch."""
        bed = tmp_path / "c.bed"
        n = callability.callable_bed(bed, CONTIG, self._calls([0, 1, 1]))
        assert n == 2
        assert bed.read_text().strip().split("\t")[1:] == ["1", "3"]

    def test_nothing_callable_writes_an_empty_bed(self, tmp_path):
        bed = tmp_path / "c.bed"
        assert callability.callable_bed(bed, CONTIG, self._calls([0, 0])) == 0
        assert bed.read_text() == ""


class TestMaskSequence:
    def _calls(self, statuses):
        return [call_position(36 if s else 0, copies=2, lambda1=18.0,
                              dispersion=float("inf")) for s in statuses]

    def test_uncallable_positions_become_n(self):
        assert callability.mask_sequence(
            "ACGTA", self._calls([1, 0, 1, 0, 1])) == "ANGNA"

    def test_all_callable_is_unchanged(self):
        assert callability.mask_sequence("ACGT", self._calls([1] * 4)) == "ACGT"

    def test_empty_sequence_is_returned_as_is(self):
        assert callability.mask_sequence("", self._calls([1])) == ""

    def test_more_calls_than_bases_does_not_overrun(self):
        """Consensus length and call count can disagree after an indel; this
        must truncate rather than raise."""
        assert callability.mask_sequence("AC", self._calls([1, 0, 0, 0])) == "AN"


class TestTrackAndSummary:
    def test_track_positions_are_one_based(self, tmp_path):
        """The track is read by humans against a genome browser; the evidence
        list is 0-based and the track must not be."""
        ev = [PositionEvidence(pos=0, depth=36, n_pass_mapq=36)]
        calls = callability.classify(ev, copies=2, lambda1=18.0,
                                     dispersion=float("inf"))
        out = tmp_path / "t.tsv"
        callability.write_track(out, CONTIG, ev, calls)
        rows = out.read_text().splitlines()
        assert rows[1].split("\t")[1] == "1"

    def test_summary_counts_reach_the_caller(self, tmp_path):
        ev = [PositionEvidence(pos=i, depth=d, n_pass_mapq=d)
              for i, d in enumerate([36, 36, 0])]
        calls = callability.classify(ev, copies=2, lambda1=18.0,
                                     dispersion=float("inf"))
        summary = callability.write_track(tmp_path / "t.tsv", CONTIG, ev, calls)
        assert summary["n_positions"] == 3
        assert summary["n_ok"] == 2
        assert summary["n_low_depth"] == 1
        # summarise() rounds to 4 dp, so compare at that precision rather than
        # asserting a float identity the function never promised.
        assert summary["callable_fraction"] == pytest.approx(2 / 3, abs=1e-4)


class TestGcByPosition:
    def test_uniform_sequence_gives_its_own_gc(self):
        assert callability.gc_by_position("GC" * 100)[50] == pytest.approx(1.0)
        assert callability.gc_by_position("AT" * 100)[50] == pytest.approx(0.0)

    def test_n_bases_are_excluded_from_the_denominator(self):
        """N is unknown, not AT. Counting it as non-GC biases the correction
        toward low GC exactly at assembly gaps."""
        assert callability.gc_by_position("N" * 50 + "G" * 50)[70] == \
            pytest.approx(1.0)

    def test_all_n_gives_zero_rather_than_dividing_by_zero(self):
        assert callability.gc_by_position("N" * 50)[25] == 0.0

    def test_one_value_per_base(self):
        assert len(callability.gc_by_position("ACGT" * 30)) == 120
