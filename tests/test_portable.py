"""Tests that `portable/lilr_cn.py` still agrees with the package.

The portable script is a deliberate duplicate: one file, samtools, no imports
from `lilrwgs`, so it can be copied onto a machine that has none of this
repository. Duplication is the price of that, and the risk it buys is drift —
a constant corrected here and not there produces a script that still runs, still
emits copy numbers, and is quietly wrong.

Nothing here needs a BAM. Every load-bearing number and every pure decision is
compared against the package's own, so a change to either side that is not made
to both fails loudly the next time pytest runs.
"""

from __future__ import annotations

import ast
import importlib.util
import sys
from pathlib import Path

import pytest

from lilrwgs import cn as cn_module
from lilrwgs import coverage as coverage_module
from lilrwgs import loci

PORTABLE = Path(__file__).resolve().parent.parent / "portable" / "lilr_cn.py"


def load_portable():
    spec = importlib.util.spec_from_file_location("lilr_cn_portable", PORTABLE)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


portable = load_portable()


class TestItIsActuallyPortable:
    def test_imports_nothing_from_this_repository(self):
        """The whole point: `scp lilr_cn.py elsewhere` has to be enough."""
        tree = ast.parse(PORTABLE.read_text())
        imported: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.update(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                imported.add(node.module or "")

        assert not any(m.startswith("lilrwgs") for m in imported), (
            "portable/lilr_cn.py imports lilrwgs, so it is not portable"
        )
        top_level = {m.split(".")[0] for m in imported if m}
        third_party = top_level - set(sys.stdlib_module_names)
        assert not third_party, (
            f"portable/lilr_cn.py imports {sorted(third_party)}; it is meant to "
            "need nothing but the standard library and samtools on PATH"
        )

    def test_no_site_specific_paths(self):
        """A path from the machine it was written on is the usual way this breaks."""
        text = PORTABLE.read_text()
        for marker in ("/wynton", "micromamba", "/scratch/", "qsub", "SGE_", "sge"):
            assert marker not in text, f"portable/lilr_cn.py mentions {marker!r}"


class TestCoordinates:
    def test_control_loci_are_the_same_loci(self):
        package = [(c.name, c.chrom, c.start, c.end, c.inside_placement)
                   for c in loci.ALL_CONTROLS]
        assert portable.CONTROLS == package

    def test_unique_windows_match(self):
        assert portable.UNIQUE_WINDOWS == loci.UNIQUE_WINDOWS

    def test_gene_spans_match(self):
        for gene, span in portable.GENE_SPANS.items():
            assert span == loci.gene_span(gene)

    def test_lilra3_intervals_and_junction_match(self):
        assert portable.LILRA3_ALT == loci.LILRA3_ALT
        assert portable.ALT_CONTIGS == loci.ALT_CONTIGS
        assert portable.LILRA3_JUNCTION == loci.LILRA3_JUNCTION
        assert (portable.LILRA3_JUNCTION_MICROHOMOLOGY
                == loci.LILRA3_JUNCTION_MICROHOMOLOGY)

    def test_slice_intervals_are_identical(self):
        """Including the merge, which is what keeps lambda_1 off double depth."""
        assert portable.slice_intervals() == loci.slice_intervals()
        assert (portable.slice_intervals(include_alts=False)
                == loci.slice_intervals(include_alts=False))

    def test_merge_absorbs_the_same_intervals(self):
        cases = [
            [("chr19", 100, 1000), ("chr19", 200, 300)],
            [("chr19", 100, 200), ("chr19", 200, 300)],
            [("chr19", 300, 400), ("chr19", 100, 200), ("chr1", 50, 60)],
        ]
        for case in cases:
            assert portable.merge_intervals(case) == loci.merge_intervals(case)


class TestThresholds:
    @pytest.mark.parametrize("name", [
        "MAPQ_STRICT", "MAPQ_ANY", "FLANK", "CHROM",
    ])
    def test_loci_constants(self, name):
        assert getattr(portable, name) == getattr(loci, name)

    @pytest.mark.parametrize("name", ["MIN_BASEQ", "MIN_Q20_LRC", "MAX_DILUTION"])
    def test_coverage_constants(self, name):
        assert getattr(portable, name) == getattr(coverage_module, name)

    @pytest.mark.parametrize("name", [
        "CN_RANGE", "AMBIGUOUS_BAND", "MIN_CLIP", "MIN_ANCHOR",
        "JUNCTION_TOLERANCE",
    ])
    def test_cn_constants(self, name):
        assert getattr(portable, name) == getattr(cn_module, name)

    def test_slice_keeps_supplementary_records(self):
        """LILRA3's evidence is supplementary; dropping it manufactures deletions."""
        from lilrwgs import extract

        assert portable.EXCLUDE_FLAGS == extract.EXCLUDE_FLAGS
        assert not portable.EXCLUDE_FLAGS & extract.SUPPLEMENTARY


class TestDecisions:
    def test_integerise_agrees_everywhere(self):
        for gene in ("LILRA3", "LILRA6", "LILRB3", "LILRB1"):
            for i in range(-100, 800):
                estimate = i / 100.0
                assert (portable.integerise(estimate, gene)
                        == cn_module.integerise(estimate, gene)), (estimate, gene)

    def test_alignment_verdict_agrees_everywhere(self):
        for q20_lrc in (0.0, 0.29, 0.3, 0.5, 0.95):
            for q20_outside in (0.0, 0.29, 0.3, 0.96):
                for dilution in (0.0, 1.0, 1.5, 1.51, 3.0):
                    assert (portable.alignment_verdict(q20_lrc, q20_outside, dilution)
                            == coverage_module.alignment_verdict(
                                q20_lrc, q20_outside, dilution))

    def test_junction_tally_agrees(self):
        chrom, left = loci.LILRA3_JUNCTION
        right = left + loci.LILRA3_JUNCTION_MICROHOMOLOGY
        lines = [
            # right-clipped at the first breakpoint end: a bearing chromosome
            f"r1\t0\t{chrom}\t{left - 49}\t60\t50M40S\t*\t0\t0\t*\t*",
            # left-clipped at the second: the same chromosome coming back out
            f"r2\t0\t{chrom}\t{right + 1}\t60\t40S50M\t*\t0\t0\t*\t*",
            # spanning: a deleted chromosome, which matches the primary assembly
            f"r3\t0\t{chrom}\t{left - 60}\t60\t120M\t*\t0\t0\t*\t*",
            # a short clip is trimming, not a breakpoint
            f"r4\t0\t{chrom}\t{left - 49}\t60\t50M4S\t*\t0\t0\t*\t*",
            "malformed",
        ]
        assert (portable.tally_junction(lines, left, right)
                == cn_module.tally_junction(lines, left, right))
        assert portable.tally_junction(lines, left, right)[0] == 2


class TestThreadBudget:
    """`--threads` is the only knob a user has, so it has to divide sensibly."""

    @pytest.mark.parametrize("threads,n_samples,expected", [
        (1, 1, (1, 1)),
        (16, 1, (1, 16)),        # one sample: the budget goes to its slice
        (16, 4, (4, 4)),
        (16, 100, (16, 1)),      # a cohort: run them side by side
        (0, 10, (1, 1)),         # nonsense in, something runnable out
    ])
    def test_allocation(self, threads, n_samples, expected):
        assert portable.allocate(threads, n_samples) == expected


class TestManifests:
    def test_reads_the_repository_manifest_format(self, tmp_path):
        manifest = tmp_path / "m.tsv"
        manifest.write_text(
            "sample_id\tcram\tcrai\tpopulation\n"
            "HG00096\thttps://example/HG00096.final.cram\thttps://x.crai\tGBR\n")
        assert portable.read_manifest(manifest) == [
            ("HG00096", "https://example/HG00096.final.cram")]

    def test_reads_a_bare_list_of_paths(self, tmp_path):
        manifest = tmp_path / "m.txt"
        manifest.write_text("/data/HG00096.slice.bam\n/data/NA12878.final.cram\n")
        assert portable.read_manifest(manifest) == [
            ("HG00096", "/data/HG00096.slice.bam"),
            ("NA12878", "/data/NA12878.final.cram")]

    def test_sample_names_lose_the_bioinformatics_suffixes(self):
        assert portable.sample_name("/a/b/HG00096.final.cram") == "HG00096"
        assert portable.sample_name("HG00096.slice.bam") == "HG00096"
        assert portable.sample_name("https://x/y/NA12878.cram?v=2") == "NA12878"


class TestRefusal:
    """A failed measurement and a true zero stay different values here too."""

    def _model(self, **kw):
        model = portable.CoverageModel("S")
        model.lambda1 = 15.0
        model.q20_lrc = 0.95
        model.alt_verdict = "alt_aware"
        for key, value in kw.items():
            setattr(model, key, value)
        return model

    def test_dead_mapq20_is_not_measured_rather_than_zero(self):
        model = self._model(q20_lrc=0.02, alt_verdict="not_alt_aware")
        call = portable._call_unique_window("S", "LILRA6", "nonexistent.bam",
                                            model, reference=None,
                                            samtools="samtools")
        assert call.status == "not_measured"
        assert call.copies is None

    def test_no_coverage_model_is_a_failure_rather_than_zero(self):
        call = portable._call_unique_window("S", "LILRB3", "nonexistent.bam",
                                            self._model(lambda1=0.0),
                                            reference=None, samtools="samtools")
        assert call.status == "failed"
        assert call.copies is None

    def test_statuses_are_the_same_three(self):
        source = PORTABLE.read_text()
        for status in ("measured", "not_measured", "failed"):
            assert f'"{status}"' in source
