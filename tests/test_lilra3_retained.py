"""LILRA3's 3' end survives its own deletion, and genotyping has to know.

The deletion removes the first 4,817 bp of the 7,126 bp calling reference. The
last 2,309 bp -- the end of intron 6, exon 7 with the stop codon, the 3' UTR and
flank -- is also on the deletion chromosome, so it sits at two copies whatever
the sample's LILRA3 copy number. Genotyping used to judge it at the gene's copy
number. In EUR50 that masked most of it as high_depth in all 12 one-copy
samples, turned the masked stop into a residue 440, and, where it slipped under
the ceiling, called it haploid over a 1:1 mix of two chromosomes -- which wrote
deletion-allele bases into hap1 in 6 of 12 samples with no flag.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

from lilrwgs import callability, genotype
from lilrwgs.callability import PositionEvidence
from lilrwgs.depth_model import Callability
from lilrwgs.loci import (CHROM, LILRA3_JUNCTION, LILRA3_RETAINED,
                          RETAINED_ANCHOR, RETAINED_ON_LOCUS_REF)

REPO = Path(__file__).resolve().parent.parent
LOCUS_REF = REPO / "resources" / "bundle" / "references" / "LILRA3_named.fa"
GFF = REPO / "resources" / "bundle" / "exon_coords" / "LILRA3.gff"
GRCH38 = (REPO / "resources" / "reference"
          / "GRCh38_full_analysis_set_plus_decoy_hla.fa")

START, END = RETAINED_ON_LOCUS_REF["LILRA3"]


def _seq(path: Path) -> str:
    return "".join(ln.strip() for ln in path.read_text().splitlines()
                   if not ln.startswith(">")).upper()


def _revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


class TestTheSegment:
    """The constant against the files it describes."""

    def test_it_is_the_tail_of_the_calling_reference(self):
        ref = _seq(LOCUS_REF)
        assert len(ref) == 7_126
        assert (START, END) == (4_817, 7_126)

    def test_it_holds_exon_7_and_the_stop_and_nothing_upstream(self):
        cds, stop = [], None
        for line in GFF.read_text().splitlines():
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            if parts[2] == "CDS":
                cds.append((int(parts[3]) - 1, int(parts[4])))
            elif parts[2] == "stop_codon":
                stop = (int(parts[3]) - 1, int(parts[4]))
        cds.sort()
        *upstream, exon7 = cds
        assert START <= exon7[0] and exon7[1] <= END
        assert stop is not None and START <= stop[0] and stop[1] <= END
        assert all(e <= START for _, e in upstream)

    def test_the_grch38_interval_ends_at_the_junction(self):
        (chrom, lo, hi), = LILRA3_RETAINED
        assert chrom == CHROM
        assert hi == LILRA3_JUNCTION[1]
        # Half-open, like everything in loci. 2,309 bp of reference against
        # 2,311 of primary: two 1 bp indels.
        assert hi - lo == 2_311

    def test_the_anchor_sits_across_the_span_start(self):
        ref = _seq(LOCUS_REF)
        anchor = RETAINED_ANCHOR["LILRA3"]
        half = len(anchor) // 2
        assert ref[START - half:START + half] == anchor

    @pytest.mark.skipif(not GRCH38.exists() or not shutil.which("samtools"),
                        reason="needs the GRCh38 analysis set and samtools")
    def test_it_is_on_the_deletion_chromosome_and_its_neighbour_is_not(self):
        """Primary chr19 *is* the deletion allele, so the retained segment is
        on it and the base before the segment is not."""
        (chrom, lo, hi), = LILRA3_RETAINED

        def primary(a, b):
            out = subprocess.run(["samtools", "faidx", str(GRCH38),
                                  f"{chrom}:{a}-{b}"],
                                 capture_output=True, text=True, check=True).stdout
            return _revcomp("".join(out.split("\n")[1:]).upper())

        ref = _seq(LOCUS_REF)
        seg = primary(lo + 1, hi)
        assert ref[START:START + 33] == seg[:33]
        assert ref[END - 100:END] == seg[-100:]
        assert ref[START - 28:START] != primary(hi + 1, hi + 28)


class TestCopiesByPosition:
    REF = _seq(LOCUS_REF)

    def test_one_copy_sample_has_two_copies_over_the_segment(self):
        out = genotype.copies_by_position("LILRA3", 1, self.REF)
        assert out[:START] == [1] * START
        assert out[START:END] == [2] * (END - START)

    def test_two_copy_sample_takes_the_old_path(self):
        assert genotype.copies_by_position("LILRA3", 2, self.REF) is None

    def test_other_genes_take_the_old_path(self):
        assert genotype.copies_by_position("LILRB1", 1, "A" * 8_803) is None

    def test_a_reference_it_was_not_measured_on_is_refused(self):
        """Shifted flanks: every coordinate still resolves, so only the anchor
        can tell. Refused at any copy number, not only where it would bite."""
        shifted = "A" * 50 + self.REF
        for copies in (1, 2):
            with pytest.raises(ValueError, match="not the one"):
                genotype.copies_by_position("LILRA3", copies, shifted)

    def test_a_truncated_reference_is_refused(self):
        with pytest.raises(ValueError, match="not the one"):
            genotype.copies_by_position("LILRA3", 1, self.REF[:5_000])


class TestClassifyAtThePresentCopyNumber:
    """Two copies' worth of reads over the segment is expected, not a pile-up."""

    LAMBDA = 15.0

    def _calls(self, copies_by_pos):
        depth = int(2 * self.LAMBDA)
        evidence = [PositionEvidence(pos=i, depth=depth, n_pass_mapq=depth)
                    for i in range(4)]
        return callability.classify(evidence, copies=1, lambda1=self.LAMBDA,
                                    dispersion=30.0, copies_by_pos=copies_by_pos)

    def test_judged_at_the_gene_copy_number_it_is_high_depth(self):
        assert {c.status for c in self._calls(None)} == {Callability.HIGH_DEPTH}

    def test_judged_at_the_copies_present_it_is_callable(self):
        calls = self._calls([1, 1, 2, 2])
        assert [c.status for c in calls] == [Callability.HIGH_DEPTH] * 2 + \
            [Callability.OK] * 2
        assert calls[2].effective_copies == 2

    def test_the_bed_can_be_split_by_ploidy(self, tmp_path):
        calls = self._calls([2, 2, 2, 2])
        n = callability.callable_bed(tmp_path / "b.bed", "LILRA3", calls,
                                     where=[False, True, True, False])
        assert n == 2
        assert (tmp_path / "b.bed").read_text() == "LILRA3\t1\t3\n"


class TestUnassignableHets:
    COPIES_BY_POS = [1] * 10 + [2] * 10     # positions 10..19 are shared

    @staticmethod
    def _rec(pos1: int, gt: str, ref: str = "A") -> str:
        return f"LILRA3\t{pos1}\t.\t{ref}\tG\t50\tPASS\t.\tGT:AD:DP\t{gt}:10,10:20"

    def test_a_het_in_the_shared_stretch_is_unassignable(self):
        body = self._rec(12, "0/1")
        assert genotype.unassignable_hets(body, 1, self.COPIES_BY_POS) == {11}

    def test_a_homozygote_there_is_carried_by_both_so_it_stands(self):
        body = self._rec(12, "1/1")
        assert genotype.unassignable_hets(body, 1, self.COPIES_BY_POS) == set()

    def test_outside_the_stretch_nothing_is_touched(self):
        body = "\n".join([self._rec(3, "1"), self._rec(5, "0/1")])
        assert genotype.unassignable_hets(body, 1, self.COPIES_BY_POS) == set()

    def test_phased_and_multi_base_records(self):
        body = "\n".join([self._rec(15, "1|0", ref="AC"), self._rec(18, "./.")])
        assert genotype.unassignable_hets(body, 1, self.COPIES_BY_POS) == {14, 15}

    def test_headers_are_ignored(self):
        body = "##fileformat=VCFv4.2\n#CHROM\tPOS\n" + self._rec(12, "0/1")
        assert genotype.unassignable_hets(body, 1, self.COPIES_BY_POS) == {11}


class TestTheStopCodonIsParsed:
    def test_parse_gff_keeps_the_stop(self):
        coords = genotype._parse_gff(GFF.read_text())
        assert coords["stop_codon"] == (5_625, 5_628)
        assert coords["exons"][-1] == (5_569, 5_628)


class TestProcessCallsTheStretchAtItsOwnPloidy:
    """genotype.process end to end, with every tool faked.

    The helpers above can each be right while process() fails to use them, so
    this pins the wiring: which HaplotypeCaller runs happen at which ploidy over
    which positions, that the parts are concatenated, and that a heterozygote in
    the stretch reaches the track, the BED and the filter as `deletion_shared`.
    """

    LAMBDA = 15.0
    HET_POS1 = 5_175        # 1-based, inside the stretch
    HOM_POS1 = 5_305

    def _run(self, tmp_path, monkeypatch, locus, copies):
        ref_fa = REPO / "resources" / "bundle" / "references" / f"{locus}_named.fa"
        ref_seq = _seq(ref_fa)
        log = {"hc": [], "concat": 0, "filter_beds": []}

        def bed_of(cmd, flag):
            return Path(cmd[cmd.index(flag) + 1]).read_text()

        class _Res:
            def __init__(self, stdout=""):
                self.stdout, self.stderr, self.returncode = stdout, "", 0

        def fake_run(cmd, **kw):
            tool = cmd[0]
            if "HaplotypeCaller" in cmd:
                log["hc"].append((int(cmd[cmd.index("-ploidy") + 1]),
                                  bed_of(cmd, "-L")))
                Path(cmd[cmd.index("-O") + 1]).touch()
            elif tool == "bcftools" and cmd[1] == "concat":
                log["concat"] += 1
                Path(cmd[cmd.index("-o") + 1]).touch()
            elif tool == "bcftools" and cmd[1] == "view" and "-T" in cmd:
                log["filter_beds"].append(bed_of(cmd, "-T"))
                return _Res("##fileformat=VCFv4.2\n")
            elif tool == "bcftools" and cmd[1] == "view" and "-o" in cmd:
                Path(cmd[cmd.index("-o") + 1]).touch()
            elif tool == "bcftools" and cmd[1] == "view" and "-H" in cmd:
                rec = "{c}\t{p}\t.\tA\tG\t50\tPASS\t.\tGT:AD:DP\t{gt}:10,10:20"
                return _Res("\n".join([
                    rec.format(c=locus, p=self.HET_POS1, gt="0/1"),
                    rec.format(c=locus, p=self.HOM_POS1, gt="1/1")]))
            elif tool == "bcftools" and cmd[1] == "consensus":
                Path(cmd[cmd.index("-o") + 1]).write_text(f">x\n{ref_seq}\n")
            elif tool == "whatshap":
                Path(cmd[cmd.index("-o") + 1]).touch()
            return _Res()

        def fake_pipeline(stages, **kw):
            sort = stages[-1]
            assert sort[:2] == ["samtools", "sort"]
            Path(sort[sort.index("-o") + 1]).write_bytes(b"BAM")
            return _Res()

        def fake_evidence(bam, n, **kw):
            per_copy = genotype.copies_by_position(locus, copies, ref_seq) \
                or [copies] * n
            return [PositionEvidence(pos=i, depth=int(self.LAMBDA * per_copy[i]),
                                     n_pass_mapq=int(self.LAMBDA * per_copy[i]))
                    for i in range(n)]

        monkeypatch.setattr(genotype, "require", lambda *a: None)
        monkeypatch.setattr(genotype, "run", fake_run)
        monkeypatch.setattr(genotype, "pipeline", fake_pipeline)
        monkeypatch.setattr(genotype.callability, "gather_evidence", fake_evidence)
        monkeypatch.setattr(genotype, "_finish_haplotype",
                            lambda s, l, hap, *a, **kw: {"hap": hap})

        reads = tmp_path / "reads"
        reads.mkdir()
        for r in ("R1", "R2"):
            (reads / f"{locus}_{r}.fq.gz").write_bytes(b"")
        result = genotype.process(
            "S", locus, reads, copies, 0, self.LAMBDA, 30.0,
            locus_ref_dir=ref_fa.parent, locus_index_dir=tmp_path,
            protein_dir=tmp_path, work_dir=tmp_path / "work",
            out_dir=tmp_path / "out")
        track = (tmp_path / "out" / "S" / f"{locus}.callability.tsv").read_text()
        status = {int(r.split("\t")[1]): r.split("\t")[-1]
                  for r in track.splitlines()[1:]}
        return result, log, status

    def test_one_copy_lilra3_is_called_at_two_ploidies_and_joined(
            self, tmp_path, monkeypatch):
        result, log, _ = self._run(tmp_path, monkeypatch, "LILRA3", 1)
        assert result.status == "ok"
        assert sorted(log["hc"]) == [(1, f"LILRA3\t0\t{START}\n"),
                                     (2, f"LILRA3\t{START}\t{END}\n")]
        assert log["concat"] == 1

    def test_a_het_in_the_stretch_is_masked_everywhere_it_counts(
            self, tmp_path, monkeypatch):
        result, log, status = self._run(tmp_path, monkeypatch, "LILRA3", 1)
        assert status[self.HET_POS1] == "deletion_shared"
        assert status[self.HOM_POS1] == "ok"
        assert result.callability["n_deletion_shared"] == 1
        # The filter ran again on a BED rebuilt without the het.
        assert len(log["filter_beds"]) == 2
        het0 = self.HET_POS1 - 1
        assert log["filter_beds"][1] == (f"LILRA3\t0\t{het0}\n"
                                         f"LILRA3\t{het0 + 1}\t{END}\n")

    def test_two_copy_lilra3_takes_the_single_call_path(self, tmp_path,
                                                        monkeypatch):
        result, log, status = self._run(tmp_path, monkeypatch, "LILRA3", 2)
        assert log["hc"] == [(2, f"LILRA3\t0\t{END}\n")]
        assert log["concat"] == 0
        assert len(log["filter_beds"]) == 1
        assert "deletion_shared" not in status.values()

    def test_another_gene_at_one_copy_takes_the_single_call_path(
            self, tmp_path, monkeypatch):
        result, log, status = self._run(tmp_path, monkeypatch, "LILRB1", 1)
        assert [p for p, _ in log["hc"]] == [1]
        assert log["concat"] == 0
        assert len(log["filter_beds"]) == 1
        assert "deletion_shared" not in status.values()
