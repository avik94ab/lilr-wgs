"""Five defects a code review found, each pinned where it was fixed.

Every one of them left output that stayed plausible, which is why none of the
existing tests saw it:

- the allele-balance filter tested the first ALT's fraction, not the minor
  allele's, so a het with three reference reads out of thirty passed;
- bowtie2 was piped into samtools by hand and its exit code never read, so a
  bowtie2 killed partway left a valid, partial BAM;
- a copy-number call that was not measured was genotyped as a measured 2, with
  nothing in the marker to say so;
- the `.alt missing` warning fired on CHM13, which correctly has no `.alt`;
- `extract` counted singletons as pairs, a bug `realign` had already fixed in
  its own copy of the parser.
"""

from __future__ import annotations

import inspect
from pathlib import Path

import pytest

from lilrwgs import extract, genotype, realign
from lilrwgs.depth_model import HET_RATIO

ROOT = Path(__file__).resolve().parent.parent


class TestAlleleBalance:
    """PING's hetRatio, on the alleles GT actually carries."""

    def test_minor_reference_allele_is_rejected(self):
        # The case the old `FMT/AD[0:1]/FMT/DP` expression passed: 27 of 30 reads
        # are ALT, so the ALT fraction is 0.9 -- but the minor allele, the
        # reference, holds 10%.
        assert not genotype.allele_balance_ok("0/1", [3, 27], 30, HET_RATIO)

    def test_minor_alt_allele_is_rejected(self):
        assert not genotype.allele_balance_ok("0/1", [27, 3], 30, HET_RATIO)

    def test_balanced_het_is_kept(self):
        assert genotype.allele_balance_ok("0|1", [15, 15], 30, HET_RATIO)

    def test_second_alt_het_is_judged_on_the_allele_it_calls(self):
        # 0/2 carries REF and ALT2. The old expression looked at ALT1, which this
        # genotype does not carry, and rejected a clean het for having 1 read.
        assert genotype.allele_balance_ok("0/2", [15, 1, 14], 30, HET_RATIO)
        assert not genotype.allele_balance_ok("0/2", [3, 1, 26], 30, HET_RATIO)

    def test_two_alt_het_checks_both_alts(self):
        assert genotype.allele_balance_ok("1/2", [0, 15, 15], 30, HET_RATIO)
        assert not genotype.allele_balance_ok("1/2", [0, 27, 3], 30, HET_RATIO)

    def test_homozygotes_and_no_calls_pass(self):
        assert genotype.allele_balance_ok("1/1", [0, 30], 30, HET_RATIO)
        assert genotype.allele_balance_ok("0/0", [30, 0], 30, HET_RATIO)
        assert genotype.allele_balance_ok("./.", [], None, HET_RATIO)

    def test_het_without_evidence_fails(self):
        # An unshown balance is exactly what the filter exists to stop.
        assert not genotype.allele_balance_ok("0/1", [], 30, HET_RATIO)
        assert not genotype.allele_balance_ok("0/1", [15, 15], None, HET_RATIO)
        assert not genotype.allele_balance_ok("0/1", [15, 15], 0, HET_RATIO)
        assert not genotype.allele_balance_ok("0/1", [15, None], 30, HET_RATIO)
        assert not genotype.allele_balance_ok("0/2", [15, 15], 30, HET_RATIO)

    def test_threshold_is_inclusive(self):
        assert genotype.allele_balance_ok("0/1", [25, 75], 100, 0.25)
        assert not genotype.allele_balance_ok("0/1", [24, 76], 100, 0.25)

    def test_record_parser_reads_fields_by_format_key(self):
        line = "c\t1\t.\tA\tC\t50\t.\t.\tGT:DP:AD:GQ\t0/1:30:3,27:99\n"
        assert not genotype._record_balanced(line, HET_RATIO)
        line = "c\t1\t.\tA\tC\t50\t.\t.\tGT:AD:DP\t0/1:14,16:30\n"
        assert genotype._record_balanced(line, HET_RATIO)
        line = "c\t1\t.\tA\tC\t50\t.\t.\tGT:AD:DP\t0/1:.:30\n"
        assert not genotype._record_balanced(line, HET_RATIO)

    def test_no_bcftools_expression_decides_balance(self):
        # The balance is decided in Python; a bcftools `-i` expression cannot
        # index AD by the alleles GT carries.
        src = inspect.getsource(genotype.process)
        assert '"-i"' not in src
        assert "_record_balanced" in src


class TestToolExitCodes:
    """bowtie2 must go through shell.pipeline, which checks every stage."""

    @pytest.mark.parametrize("func", [genotype.process])
    def test_genotype_does_not_hand_pipe_bowtie2(self, func):
        src = inspect.getsource(func)
        assert "subprocess.Popen" not in src
        assert "pipeline(" in src

    def test_process_sample_does_not_hand_pipe_bowtie2(self):
        src = (ROOT / "scripts" / "process_sample.py").read_text()
        assert "subprocess.Popen" not in src
        assert "stderr=subprocess.DEVNULL" not in src


class TestCnStatus:
    """A fallback 2 is not a measured 2, and the marker says which."""

    @pytest.fixture
    def cn_tsv(self, tmp_path):
        p = tmp_path / "cn.tsv"
        p.write_text(
            "sample\tgene\tcopies\tstatus\n"
            "S1\tLILRA3\t1\tmeasured\n"
            "S1\tLILRA6\t\tnot_measured\n"
            "S1\tLILRB3\t3\tmeasured\n"
            "S2\tLILRA3\t\tfailed\n"
            "S2\tLILRA6\t\tmeasured\n"
        )
        return str(p)

    def test_measured(self, cn_tsv):
        assert genotype._copies_from(cn_tsv, "S1", "LILRA3") == (1, 0, "measured")

    def test_not_measured_falls_back_to_two_and_says_so(self, cn_tsv):
        copies, partner, status = genotype._copies_from(cn_tsv, "S1", "LILRA6")
        assert (copies, status) == (2, "not_measured")
        assert partner == 3

    def test_failed(self, cn_tsv):
        assert genotype._copies_from(cn_tsv, "S2", "LILRA3")[::2] == (2, "failed")

    def test_measured_without_a_number_is_failed(self, cn_tsv):
        assert genotype._copies_from(cn_tsv, "S2", "LILRA6")[::2] == (2, "failed")

    def test_gene_without_a_call_is_not_called(self, cn_tsv):
        assert genotype._copies_from(cn_tsv, "S1", "LILRB1")[::2] == (2, "not_called")

    def test_status_reaches_the_marker_and_the_summary(self):
        r = genotype.GenotypeResult(sample="S", locus="LILRA3", copies=2,
                                    cn_status="failed")
        assert r.as_dict()["cn_status"] == "failed"
        src = (ROOT / "scripts" / "summarise.py").read_text()
        assert '"cn_status"' in src


class TestAltWarning:
    """Only a target with ALT contigs needs the `.alt` file."""

    def _stats_for(self, monkeypatch, tmp_path, **kwargs):
        seen = []

        class Recording(realign.RealignStats):
            def __init__(self, *a, **k):
                super().__init__(*a, **k)
                seen.append(self)

        class Stop(Exception):
            pass

        def stop(*_a, **_k):
            raise Stop

        monkeypatch.setattr(realign, "RealignStats", Recording)
        monkeypatch.setattr(realign, "index_is_alt_aware", lambda _i: False)
        monkeypatch.setattr(realign, "read_contigs", stop)
        with pytest.raises(Stop):
            realign.realign_sample("S", "in.bam", "ref.fa", tmp_path / "idx",
                                   tmp_path / "out.bam", tmpdir=tmp_path, **kwargs)
        return seen[0]

    def test_grch38_target_warns(self, monkeypatch, tmp_path):
        stats = self._stats_for(monkeypatch, tmp_path)
        assert any(".alt is missing" in w for w in stats.warnings)

    def test_chm13_target_does_not(self, monkeypatch, tmp_path):
        stats = self._stats_for(monkeypatch, tmp_path, expect_alt=False)
        assert not any(".alt" in w for w in stats.warnings)
        assert stats.alt_aware is False      # still recorded, just not alarmed

    def test_lilra6_cn_passes_the_target_assembly(self):
        src = (ROOT / "scripts" / "lilra6_cn.py").read_text()
        assert 'expect_alt=target_assembly == "GRCh38"' in src


class TestFastqCountsInExtract:
    STDERR = ("[M::bam2fq_mainloop] discarded 2362 singletons\n"
              "[M::bam2fq_mainloop] processed 167082 reads\n")

    def test_singletons_are_not_counted_as_pairs(self):
        assert extract._parse_fastq_counts(self.STDERR) == (82360, 2362)

    def test_one_parser_serves_both_front_ends(self):
        assert realign._parse_fastq_counts is extract._parse_fastq_counts
