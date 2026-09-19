"""Tests for CDS extraction and IUPAC-aware translation.

`sequences.py` is vendored from `lilr-genotyper` and only partly live: this
pipeline reaches `extract_sequences`, `_translate` and `_revcomp`, while
`build_consensus` and its `MIN_CONS_DEPTH` belong to the upstream pileup path
and are not used. That boundary is load-bearing rather than incidental — masking
here is driven by `depth_model`, and a fixed depth constant reappearing in the
consensus would be the predecessor's failure mode returning by the back door
(`docs/variant_calling.md` §3.4). `TestTheDeadCodeStaysDead` pins it.

Seven of the eleven LILR genes are on the minus strand, so the strand handling
in `extract_sequences` is the common path, not the exception. Getting it wrong
produces a protein — just the wrong one.
"""

from __future__ import annotations

import pytest

from lilrwgs import sequences
from lilrwgs.sequences import _revcomp, _translate, extract_sequences


class TestRevcomp:
    def test_plain_bases(self):
        assert _revcomp("ACGT") == "ACGT"
        assert _revcomp("AAAA") == "TTTT"

    def test_case_is_preserved(self):
        assert _revcomp("acgt") == "acgt"

    @pytest.mark.parametrize("base,complement", [
        ("M", "K"), ("K", "M"),      # A/C  <-> T/G
        ("R", "Y"), ("Y", "R"),      # A/G  <-> T/C
        ("W", "W"), ("S", "S"),      # A/T and C/G are self-complementary
        ("V", "B"), ("B", "V"),
        ("H", "D"), ("D", "H"),
        ("N", "N"),
    ])
    def test_iupac_ambiguity_codes_complement_correctly(self, base, complement):
        """A heterozygous position is written as an IUPAC code, and seven of the
        eleven genes are minus-strand — so a wrong complement here silently
        swaps which allele a heterozygote is reported to carry."""
        assert _revcomp(base) == complement

    def test_round_trip_is_identity(self):
        seq = "ACGTRYKMWSBDHVN"
        assert _revcomp(_revcomp(seq)) == seq


class TestTranslate:
    def test_a_simple_orf(self):
        # ATG AAA GGG -> M K G
        assert _translate("ATGAAAGGG") == "MKG"

    def test_translation_stops_at_a_stop_codon(self):
        assert _translate("ATGAAATAAGGG") == "MK"

    def test_a_trailing_partial_codon_is_ignored(self):
        """Not an error: a consensus can end mid-codon after an indel."""
        assert _translate("ATGAAAGG") == "MK"

    def test_n_gives_x_rather_than_a_guess(self):
        assert _translate("ATGNNNGGG") == "MXG"

    def test_ambiguity_that_is_synonymous_resolves(self):
        """CTN is leucine whatever N is, so the ambiguity does not reach the
        protein. Emitting X here would discard real information."""
        assert _translate("CTN") == "L"

    def test_ambiguity_that_is_not_synonymous_gives_x(self):
        """ARG spans AAG (K) and AGG (R); the honest answer is X."""
        assert _translate("ARG") == "X"

    def test_an_ambiguous_codon_that_is_always_stop_still_stops(self):
        # TAR covers TAA and TAG, both stop.
        assert _translate("ATGTAR") == "M"

    def test_empty_input(self):
        assert _translate("") == ""


class TestExtractSequences:
    """`coords` uses 0-based half-open exons sorted low->high, and `strand`
    decides both the order they are joined in and whether each is complemented."""

    def _plus(self):
        # consensus:  0.........10........20
        #             ATGAAAGGGTTTCCCTAAGGG
        return {
            "mrna_start": 0, "mrna_end": 18, "strand": "+",
            "exons": [(0, 9), (9, 18)],
        }

    def test_plus_strand_joins_exons_in_coordinate_order(self):
        cons = "ATGAAAGGG" "TTTCCCTAA" "GGG"
        gdna, cdna, protein = extract_sequences(cons, self._plus())
        assert cdna == "ATGAAAGGGTTTCCCTAA"
        assert protein == "MKGFP"          # stops at TAA

    def test_plus_strand_gdna_is_the_mrna_span(self):
        cons = "ATGAAAGGG" "TTTCCCTAA" "GGG"
        gdna, _, _ = extract_sequences(cons, self._plus())
        assert gdna == "ATGAAAGGGTTTCCCTAA"

    def test_minus_strand_reverses_exon_order_and_complements_each(self):
        """The common case: seven of eleven LILR genes are minus-strand.

        Built by reverse-complementing a known plus-strand transcript, so the
        expected protein is known independently of the function under test.
        """
        plus_cdna = "ATGAAAGGGTTTCCCTAA"
        cons = _revcomp(plus_cdna)          # 18 bp, minus-strand consensus
        coords = {
            "mrna_start": 0, "mrna_end": 18, "strand": "-",
            "exons": [(0, 9), (9, 18)],
        }
        gdna, cdna, protein = extract_sequences(cons, coords)
        assert cdna == plus_cdna
        assert protein == "MKGFP"

    def test_minus_strand_gdna_is_reverse_complemented(self):
        plus = "ATGAAAGGGTTTCCCTAA"
        cons = _revcomp(plus)
        coords = {"mrna_start": 0, "mrna_end": 18, "strand": "-",
                  "exons": [(0, 18)]}
        gdna, _, _ = extract_sequences(cons, coords)
        assert gdna == plus

    def test_introns_are_excluded_from_cdna_but_not_from_gdna(self):
        """The distinction the two outputs exist to make."""
        cons = "ATGAAA" "TTTTTT" "GGGTAA"      # exon, intron, exon
        coords = {"mrna_start": 0, "mrna_end": 18, "strand": "+",
                  "exons": [(0, 6), (12, 18)]}
        gdna, cdna, protein = extract_sequences(cons, coords)
        assert cdna == "ATGAAAGGGTAA"
        assert "TTTTTT" in gdna
        assert protein == "MKG"

    def test_strand_defaults_to_plus_when_absent(self):
        coords = {"mrna_start": 0, "mrna_end": 9, "exons": [(0, 9)]}
        _, cdna, _ = extract_sequences("ATGAAAGGG", coords)
        assert cdna == "ATGAAAGGG"

    def test_masked_positions_propagate_to_the_protein_as_x(self):
        """An N from the callability mask must not become a silent amino acid.

        This is the end of the chain that `mask_sequence` starts: a position the
        depth model refused becomes N in the consensus, X in the protein, and is
        never quietly resolved to whatever the reference happened to say.
        """
        cons = "ATG" + "NNN" + "GGGTAA"
        coords = {"mrna_start": 0, "mrna_end": 12, "strand": "+",
                  "exons": [(0, 12)]}
        _, _, protein = extract_sequences(cons, coords)
        assert protein == "MXG"


class TestTheDeadCodeStaysDead:
    """`build_consensus` masks at a fixed `MIN_CONS_DEPTH` and flags
    low-confidence at another fixed constant. Both belong to the predecessor's
    pileup path. Wiring either into this pipeline would put a second, constant
    masking threshold beside the callable track — which is the exact pattern
    `docs/variant_calling.md` §3.4 rules out, and the reason a position with
    depth in [15, 20) used to survive the mask but lose its call.
    """

    def test_the_live_surface_is_what_genotype_imports(self):
        import inspect

        from lilrwgs import genotype
        src = inspect.getsource(genotype)
        assert "from .sequences import extract_sequences" in src

    def test_genotype_does_not_reach_the_constant_masking_path(self):
        import inspect

        from lilrwgs import genotype
        src = inspect.getsource(genotype)
        for forbidden in ("build_consensus", "MIN_CONS_DEPTH",
                          "LOW_CONF_DEPTH", "find_novel_snps"):
            assert forbidden not in src, (
                f"{forbidden} reached the genotyping path; masking here is "
                "driven by depth_model, and a fixed depth constant in the "
                "consensus is the predecessor's failure mode returning")

    def test_the_constants_still_exist_and_are_still_fixed(self):
        """Not a recommendation to use them — a record that they are constants,
        so that if anyone does wire them up the diff shows what it costs."""
        assert isinstance(sequences.MIN_CONS_DEPTH, int)
        assert isinstance(sequences.LOW_CONF_DEPTH, int)
