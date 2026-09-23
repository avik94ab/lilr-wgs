"""The gene model comes from the unmasked consensus; the sequence from the masked one.

`genotype._finish_haplotype` finds CDS coordinates with miniprot. It used to run
miniprot on the *masked* consensus, and an N-run over a splice site reshapes the
model: on HG00097 a masked intron boundary in LILRA3 merged exons 3 and 4 and
translated 147 bp of intron (488 aa against a true 439), and a masked final exon
in LILRB3 dropped 64 residues instead of reporting them as X. A masked codon is
an unknown residue, not a different gene structure.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

from lilrwgs import genotype

REPO = Path(__file__).resolve().parent.parent
LOCUS_REFS = REPO / "resources" / "bundle" / "references"
PROTEINS = REPO / "resources" / "bundle" / "exon_coords"

# The N-run over LILRA3's intron-3/exon-4 boundary in HG00097's hap1, 1-based.
HG00097_LILRA3_NRUN = (2307, 2397)


def _mask(seq: str, start: int, end: int) -> str:
    """Mask 1-based inclusive [start, end]."""
    return seq[:start - 1] + "N" * (end - start + 1) + seq[end:]


class TestAnnotatesUnmasked:
    """Pure: miniprot is faked, so this pins which sequence it is shown."""

    # A 2-exon plus-strand gene: ATG AAA | intron | TTT GGG TAA
    UNMASKED = "CC" + "ATGAAA" + "GTAAGTTTTTTTTTTTTTCAG" + "TTTGGGTAA" + "CC"
    GFF = ("x\tminiprot\tmRNA\t3\t38\t.\t+\t.\t.\n"
           "x\tminiprot\tCDS\t3\t8\t.\t+\t.\t.\n"
           "x\tminiprot\tCDS\t30\t38\t.\t+\t.\t.\n")

    def _run(self, tmp_path, monkeypatch, masked, unmasked):
        shown = []
        real_run = subprocess.run

        def fake_run(cmd, *a, **kw):
            if cmd and cmd[0] == "miniprot":
                shown.append(Path(cmd[2]).read_text().split("\n")[1])
                return subprocess.CompletedProcess(cmd, 0, stdout=self.GFF, stderr="")
            return real_run(cmd, *a, **kw)

        monkeypatch.setattr(genotype.subprocess, "run", fake_run)
        out = tmp_path / "out"
        out.mkdir()
        row = genotype._finish_haplotype("S", "G", 1, masked, tmp_path, out,
                                         tmp_path, "miniprot", annotate_on=unmasked)
        return shown, row, out

    def test_miniprot_sees_the_unmasked_consensus(self, tmp_path, monkeypatch):
        # Mask the acceptor splice site and the first codon of exon 2.
        masked = _mask(self.UNMASKED, 27, 32)
        shown, _, _ = self._run(tmp_path, monkeypatch, masked, self.UNMASKED)
        assert shown == [self.UNMASKED]

    def test_masked_codon_is_x_not_a_different_length(self, tmp_path, monkeypatch):
        masked = _mask(self.UNMASKED, 27, 32)
        _, row, out = self._run(tmp_path, monkeypatch, masked, self.UNMASKED)
        protein = (out / "G_hap1_protein.fasta").read_text().split("\n")[1]
        assert protein == "MKXG"
        assert row["protein_len"] == 4
        # The written gDNA is still the masked one.
        gdna = (out / "G_hap1_gdna.fasta").read_text().split("\n")[1]
        assert gdna == masked

    def test_length_mismatch_is_refused(self, tmp_path, monkeypatch):
        with pytest.raises(ValueError, match="differ in length"):
            self._run(tmp_path, monkeypatch, self.UNMASKED[:-1], self.UNMASKED)


@pytest.mark.skipif(shutil.which("miniprot") is None, reason="miniprot not on PATH")
class TestLILRA3SpliceSiteMask:
    """Real miniprot, shipped LILRA3 reference, HG00097's actual N-run."""

    def _protein(self, tmp_path, masked, unmasked):
        out = tmp_path / "out"
        out.mkdir(exist_ok=True)
        genotype._finish_haplotype("HG00097", "LILRA3", 1, masked, tmp_path, out,
                                   PROTEINS, "miniprot", annotate_on=unmasked)
        return (out / "LILRA3_hap1_protein.fasta").read_text().split("\n")[1]

    def test_masked_splice_site_keeps_the_true_length(self, tmp_path):
        ref = genotype._read_fasta(LOCUS_REFS / "LILRA3_named.fa").upper()
        masked = _mask(ref, *HG00097_LILRA3_NRUN)
        truth = self._protein(tmp_path, ref, ref)
        called = self._protein(tmp_path, masked, ref)
        assert len(truth) == len(called) == 439
        assert "X" in called
        # Everywhere it is resolved, the called protein is the true one.
        assert all(c == t for c, t in zip(called, truth) if c != "X")

    def test_annotating_the_masked_sequence_is_what_broke(self, tmp_path):
        """The old behaviour, kept as evidence that the test above can fail."""
        ref = genotype._read_fasta(LOCUS_REFS / "LILRA3_named.fa").upper()
        masked = _mask(ref, *HG00097_LILRA3_NRUN)
        assert len(self._protein(tmp_path, masked, masked)) != 439
