"""Tests for the extract-and-realign front end.

The assembly tests carry the weight. Chromosome 19 is called `chr19` in both
GRCh38 and GRCh37, so a GRCh37 CRAM passes every contig-name check this package
had before `source_assembly` existed — and then returns a different
half-megabase of chromosome 19, realigns it perfectly, and reports copy numbers
for whatever genes happen to live there. Nothing downstream can tell that from a
real answer, which is why the refusal is tested and not just written.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

from lilrwgs import loci, realign  # noqa: E402
from lilrwgs.shell import ToolError  # noqa: E402

GRCH38 = {"chr19": 58_617_616, "chr1": 248_956_422}
GRCH37 = {"chr19": 59_128_983, "chr1": 249_250_621}
GRCH38_ALTS = {**GRCH38, **{c: 1_000_000 for c in loci.ALT_CONTIGS}}


class TestSourceAssembly:
    def test_grch38_is_accepted(self):
        assert realign.source_assembly(GRCH38) == ("GRCh38", "chr19")

    def test_ensembl_naming_is_resolved_not_refused(self):
        """Same assembly, spelled differently; the coordinates still apply."""
        assembly, name = realign.source_assembly({"19": 58_617_616})
        assert (assembly, name) == ("GRCh38", "19")

    def test_grch37_is_refused_by_length(self):
        with pytest.raises(ToolError) as exc:
            realign.source_assembly(GRCH37)
        assert "GRCh37" in str(exc.value)

    def test_grch37_refusal_explains_the_name_collision(self):
        """The message has to say why a name check was not enough."""
        with pytest.raises(ToolError) as exc:
            realign.source_assembly(GRCH37)
        assert "same name" in str(exc.value)

    def test_unknown_length_is_refused_rather_than_assumed(self):
        with pytest.raises(ToolError) as exc:
            realign.source_assembly({"chr19": 12_345})
        assert "matches no assembly" in str(exc.value)

    def test_missing_chr19_is_refused(self):
        with pytest.raises(ToolError) as exc:
            realign.source_assembly({"chr1": 248_956_422})
        assert "chromosome 19" in str(exc.value)


class TestResolveRegions:
    def test_alt_contigs_included_when_present(self):
        regions, has_alts = realign.resolve_regions(GRCH38_ALTS, "chr19")
        assert has_alts
        assert any("_alt" in r for r in regions)

    def test_alt_contigs_dropped_when_absent(self):
        """A primary-only reference costs LILRA3's depth route and nothing else."""
        regions, has_alts = realign.resolve_regions(GRCH38, "chr19")
        assert not has_alts
        assert not any("_alt" in r for r in regions)
        assert regions

    def test_contig_name_follows_the_header(self):
        regions, _ = realign.resolve_regions({"19": 58_617_616}, "19")
        assert regions
        assert all(r.startswith("19:") for r in regions)

    def test_coordinates_do_not_move_with_the_name(self):
        """Only the spelling is the input's; the intervals are loci.py's."""
        ucsc, _ = realign.resolve_regions(GRCH38, "chr19")
        ensembl, _ = realign.resolve_regions({"19": 58_617_616}, "19")
        assert [r.split(":", 1)[1] for r in ucsc] == \
               [r.split(":", 1)[1] for r in ensembl]

    def test_intervals_are_clipped_to_the_contig(self):
        """A flanked interval can run past the end, which samtools takes as an
        error rather than as a truncation."""
        short = {"chr19": 54_300_000}
        regions, _ = realign.resolve_regions(short, "chr19")
        for region in regions:
            end = int(region.rsplit("-", 1)[1])
            assert end <= 54_300_000

    def test_every_control_locus_survives(self):
        """lambda_1 comes from the controls; a region set without them realigns
        beautifully and has no baseline to divide by."""
        regions, _ = realign.resolve_regions(GRCH38_ALTS, "chr19")
        spans = [(int(r.split(":")[1].split("-")[0]), int(r.rsplit("-", 1)[1]))
                 for r in regions if r.startswith("chr19:")]
        for control in loci.ALL_CONTROLS:
            assert any(s <= control.start + 1 and e >= control.end
                       for s, e in spans), f"{control.name} is not covered"


class TestAltAwareness:
    def test_missing_alt_file_is_detected(self, tmp_path):
        assert not realign.index_is_alt_aware(tmp_path / "ref.fa")

    def test_empty_alt_file_counts_as_missing(self, tmp_path):
        (tmp_path / "ref.fa.alt").write_text("")
        assert not realign.index_is_alt_aware(tmp_path / "ref.fa")

    def test_present_alt_file_is_detected(self, tmp_path):
        (tmp_path / "ref.fa.alt").write_text("chr19_GL949746v1_alt\t0\n")
        assert realign.index_is_alt_aware(tmp_path / "ref.fa")


class TestScratchCwd:
    """bwa and samtools run with cwd set to a scratch dir, so every path handed
    to them has to be absolute first.

    A relative index base cost a whole submitted job: bwa reports `fail locate
    index files`, which names the index and says nothing about the working
    directory, so it reads as a broken or incomplete index rather than as a
    resolution problem.
    """

    def test_align_resolves_relative_paths(self, tmp_path, monkeypatch):
        seen: list[str] = []

        def fake_pipeline(stages, **kwargs):
            seen.extend(stages[0])
            raise RuntimeError("stop here; the command line is what is tested")

        monkeypatch.setattr(realign, "pipeline", fake_pipeline)
        monkeypatch.setattr(realign, "require", lambda *a: None)
        monkeypatch.chdir(tmp_path)

        with pytest.raises(RuntimeError):
            realign.align("S", "r1.fq.gz", "r2.fq.gz", "ref.fa",
                          tmp_path / "out.bam", tmpdir=str(tmp_path))

        relative = [a for a in seen
                    if a.endswith(("ref.fa", "r1.fq.gz", "r2.fq.gz"))
                    and not Path(a).is_absolute()]
        assert not relative, f"left relative: {relative}"


class TestBwaArgs:
    def test_soft_clipping_is_requested(self):
        """-Y is load-bearing: LILRA3's evidence arrives as supplementary
        records on the alt contigs, and hard-clipped they carry no sequence."""
        assert "-Y" in realign.BWA_ARGS

    def test_batch_size_is_fixed(self):
        """Without -K the output depends on the thread count."""
        assert "-K" in realign.BWA_ARGS


class TestLilra3RouteOnRealignedBams:
    """A realigned BAM cannot support LILRA3's alt-contig depth route.

    That route counts MAPQ-0 depth *including supplementary records*, and on
    HG00138's CRAM slice 135 of the 964 reads with an alt-contig record have no
    primary record in the LRC at all — their primaries are on chr2, chr3, chrX
    and elsewhere, repeat-derived reads with a supplementary hit on the LILRA3
    contigs. An extraction of the LRC cannot contain them, and the FASTQ step
    drops supplementary records because emitting one would write a read twice.

    Measured cost: 1,825 alt-contig records as-is against 1,448 realigned on
    HG00138 (-21%), and across 88 samples a median 0.404 copies low, never high,
    converting 19 true CN 2 calls to CN 1. The junction assay is unaffected and
    is what the realign path uses instead.
    """

    def test_alt_depth_is_skipped_when_the_bam_is_a_realignment(self, monkeypatch):
        from lilrwgs import cn

        called: list = []
        monkeypatch.setattr(cn, "_mean_depth",
                            lambda *a, **k: called.append(a) or (10.0, 100))
        monkeypatch.setattr(cn, "junction_counts", lambda *a, **k: (40, 0))

        model = _model_with_controls()
        call = cn._call_lilra3("S", "x.bam", model, reference=None,
                               samtools="samtools", alt_depth_valid=False)

        assert not called, "the alt-contig depth route ran on a realigned BAM"
        assert call.method == "junction"
        assert call.estimate == pytest.approx(2.0)

    def test_alt_depth_is_used_on_a_cram_slice(self, monkeypatch):
        """The default path must not change: a slice straight out of a CRAM does
        hold those reads, and the depth route is the validated one there."""
        from lilrwgs import cn

        monkeypatch.setattr(cn, "_mean_depth", lambda *a, **k: (10.0, 100))
        monkeypatch.setattr(cn, "junction_counts", lambda *a, **k: (40, 0))

        call = cn._call_lilra3("S", "x.bam", _model_with_controls(),
                               reference=None, samtools="samtools")
        assert call.method == "alt_depth_q0"

    def test_the_note_says_which_reason_applies(self, monkeypatch):
        """Two different reasons the same route is unavailable; a reader of the
        output has to be able to tell them apart."""
        from lilrwgs import cn

        monkeypatch.setattr(cn, "junction_counts", lambda *a, **k: (40, 0))
        call = cn._call_lilra3("S", "x.bam", _model_with_controls(),
                               reference=None, samtools="samtools",
                               alt_depth_valid=False)
        assert "regional extraction" in " ".join(call.notes)


def _model_with_controls():
    """A CoverageModel with enough MAPQ-0 control depth to normalise on."""
    from lilrwgs.coverage import ControlMeasurement, CoverageModel

    model = CoverageModel(sample="S")
    model.lambda1 = 15.0
    model.controls = [
        ControlMeasurement(name=c.name, chrom=c.chrom, start=c.start, end=c.end,
                           inside_placement=c.inside_placement,
                           mean_q0=30.0, mean_q20=30.0, n_bases=1000)
        for c in loci.CONTROL_LOCI
    ]
    return model


class TestFastqCounts:
    """`processed` includes singletons, so `processed / 2` is not the pair count.

    Measured on HG00119: 167,082 processed, 2,362 singletons, and 82,360 records
    in each of R1 and R2. The error lands in the denominator of the singleton
    rate that decides whether to warn, which is the only reason it matters.
    """

    STDERR = ("[M::bam2fq_mainloop] discarded 2362 singletons\n"
              "[M::bam2fq_mainloop] processed 167082 reads\n")

    def test_singletons_are_not_counted_as_pairs(self):
        pairs, singletons = realign._parse_fastq_counts(self.STDERR)
        assert (pairs, singletons) == (82_360, 2_362)

    def test_unparseable_stderr_returns_zeros(self):
        """A miss must not look like a successful run of zero reads; the caller
        treats zero pairs as an error on its own terms."""
        assert realign._parse_fastq_counts("something else entirely") == (0, 0)


class TestInputList:
    def _read(self, tmp_path, text):
        from lilra6_cn import read_inputs
        p = tmp_path / "inputs.tsv"
        p.write_text(text)
        return read_inputs(p)

    def test_one_column_derives_the_sample_name(self, tmp_path):
        rows = self._read(tmp_path, "/data/HG00096.final.cram\n")
        assert rows == [{"sample": "HG00096", "source": "/data/HG00096.final.cram",
                         "index": None}]

    def test_two_columns_take_the_index(self, tmp_path):
        rows = self._read(tmp_path, "/d/HG1.cram\t/i/HG1.crai\n")
        assert rows[0]["index"] == "/i/HG1.crai"

    def test_three_columns_take_an_explicit_name(self, tmp_path):
        rows = self._read(tmp_path, "NA12878\t/d/odd.cram\t/i/odd.crai\n")
        assert rows[0]["sample"] == "NA12878"

    def test_urls_keep_their_sample_name(self, tmp_path):
        rows = self._read(tmp_path, "https://x/y/HG00097.final.cram\n")
        assert rows[0]["sample"] == "HG00097"

    def test_staged_slices_keep_their_sample_name(self, tmp_path):
        rows = self._read(tmp_path, "/scratch/HG00119.slice.bam\n")
        assert rows[0]["sample"] == "HG00119"

    def test_comments_and_blanks_are_skipped(self, tmp_path):
        rows = self._read(tmp_path, "# a note\n\n/d/HG1.cram\n")
        assert len(rows) == 1

    def test_header_line_is_skipped(self, tmp_path):
        rows = self._read(tmp_path, "sample_id\tcram\nHG1\t/d/x.cram\n")
        assert len(rows) == 1 and rows[0]["sample"] == "HG1"

    def test_duplicate_sample_names_are_refused(self, tmp_path):
        """Two rows with one name means one silently overwrites the other."""
        with pytest.raises(SystemExit):
            self._read(tmp_path, "/a/HG1.cram\n/b/HG1.cram\n")

    def test_empty_list_is_refused(self, tmp_path):
        with pytest.raises(SystemExit):
            self._read(tmp_path, "# nothing here\n")
