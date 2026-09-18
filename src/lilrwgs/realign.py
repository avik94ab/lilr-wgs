"""Extract the LRC, then realign it — making the alignment ours, not the input's.

Every MAPQ-20 number this pipeline reads is conditional on the input CRAM having
been aligned to GRCh38 *with the ``.alt`` file*. For the 1000 Genomes 30x set
that holds, because NYGC ran `bwa mem -Y` against the analysis set with its alt
index. For a CRAM from anywhere else it is an assumption, and
:mod:`lilrwgs.coverage` can only detect that it was violated — turning LILRA6
and LILRB3 into `not_measured` — never repair it.

This module repairs it. The reads over the LRC and the control loci are pulled
out of whatever alignment they arrived in, converted back to FASTQ, and realigned
against the GRCh38 analysis set with its alt index. Downstream, nothing changes:
:func:`lilrwgs.coverage.measure` and :func:`lilrwgs.cn.call_sample` take a BAM in
GRCh38 coordinates and do not care which aligner produced it. What changes is
that ALT-awareness becomes a property of *this* pipeline, asserted once at
install time by the presence of the ``.alt`` file, rather than a property of
every input file that has to be measured per sample and can only be refused.

**The extraction is the whole slice, not the LILR cluster.** λ₁ — the depth one
haploid copy yields, which every threshold and every copy number is expressed in
— is measured at control loci a few hundred kilobases outside the cluster. An
extraction shaped for the LILR genes alone would realign beautifully and have no
baseline to divide by. :func:`lilrwgs.loci.slice_intervals` already collects the
union for exactly this reason, and it is used here unchanged.

**Realignment is genome-wide, against the whole reference.** Aligning LRC reads
against an LRC-sized index would be faster and would invent MAPQ: a read whose
true home is a decoy or a KIR contig has nowhere else to go in a small index, so
it lands in the LRC at MAPQ 60 and inflates the very windows being measured.
The index is large and the read set is small, so the cost is the index load, not
the alignment.

**What this does not fix.** Duplicates are still the input's judgement — they are
excluded on the flags that arrive and never recomputed, because duplicate
detection needs the whole library and a 550 kb slice is not it. And a read the
input aligner put somewhere entirely outside these intervals is not here to be
rescued; extraction can only be as complete as the alignment it reads from,
which is why :func:`source_assembly` refuses an assembly whose coordinates
these intervals do not describe rather than extracting the wrong half-megabase.
"""

from __future__ import annotations

import os
import shutil
import time
from dataclasses import dataclass, field
from pathlib import Path

from . import loci, loci_chm13
from .extract import SUPPLEMENTARY, cram_env, read_contigs
from .shell import ToolError, pipeline, require, run

# Assemblies told apart by the length of chromosome 19, which they share the
# name of and differ on by half a megabase. Checking a length rather than the
# ``@SQ AS:`` tag is deliberate: AS is optional, frequently absent, and
# frequently wrong, whereas LN is structural and has to be right for the file to
# be readable at all.
CHR19_LENGTH = {
    58_617_616: "GRCh38",
    59_128_983: "GRCh37",
    61_707_364: "CHM13v2.0",
}

# Which coordinate table describes each assembly. An assembly this pipeline can
# read but has no table for is refused rather than approximated -- GRCh37 is in
# CHR19_LENGTH so that it can be *named* in the refusal, not so that it can be
# used.
LOCI_BY_ASSEMBLY = {
    "GRCh38": loci,
    "CHM13v2.0": loci_chm13,
}

# Contig-name spellings for chromosome 19, in the order they are tried. UCSC
# style is what the analysis set and this package use; Ensembl style is the same
# assembly spelled differently and is common enough in CRAMs from outside the
# 1000 Genomes tree to be worth resolving rather than refusing.
CHR19_ALIASES = ["chr19", "19"]

# bwa settings. Matched to what NYGC ran for the 1000 Genomes 30x set, because
# the point of realigning against EBI's index is to reproduce their placements
# rather than to resemble them:
#
#   -Y  soft-clip supplementary alignments. Load-bearing, not cosmetic: LILRA3
#       is absent from the primary assembly, so on the four alt contigs that
#       carry it essentially all of its evidence arrives as supplementary
#       records, and `lilrwgs.cn` counts them. Hard-clipped (the default) they
#       still exist but carry no sequence, and the depth route reads short.
#   -K  fix the batch size, so the output does not depend on the thread count.
#       Without it, two runs of the same sample at different --threads produce
#       different insert-size estimates and, occasionally, different MAPQ.
BWA_ARGS = ["-Y", "-K", "100000000"]


@dataclass
class RealignStats:
    """What was realigned, and enough context to know whether to believe it."""

    sample: str
    source: str
    source_assembly: str = ""
    chr19_name: str = ""
    n_source_records: int = 0
    n_pairs: int = 0
    n_singletons: int = 0
    n_realigned: int = 0
    had_alt_contigs: bool = False
    alt_aware: bool = False
    slice_s: float = 0.0
    align_s: float = 0.0
    warnings: list[str] = field(default_factory=list)

    def as_row(self) -> dict:
        return {
            "sample": self.sample,
            "source": self.source,
            "source_assembly": self.source_assembly,
            "chr19_name": self.chr19_name,
            "n_source_records": self.n_source_records,
            "n_pairs": self.n_pairs,
            "n_singletons": self.n_singletons,
            "n_realigned": self.n_realigned,
            "had_alt_contigs": self.had_alt_contigs,
            "alt_aware": self.alt_aware,
            "slice_s": round(self.slice_s, 1),
            "align_s": round(self.align_s, 1),
            "warnings": ";".join(self.warnings),
        }


def source_assembly(contigs: dict[str, int]) -> tuple[str, str]:
    """``(assembly, chr19_name)`` for a header, or raise saying why not.

    Refusing is the point. Every interval in :mod:`lilrwgs.loci` is a GRCh38
    coordinate, and chromosome 19 exists under the same name in GRCh37 — so an
    unchecked extraction from a GRCh37 CRAM does not fail, it quietly returns a
    different half-megabase of chromosome 19, realigns it, and reports copy
    numbers for whatever happened to be there. That is the failure mode this
    package treats as the most dangerous available: output that stays plausible.
    """
    for name in CHR19_ALIASES:
        if name in contigs:
            assembly = CHR19_LENGTH.get(contigs[name], "")
            if assembly in LOCI_BY_ASSEMBLY:
                return assembly, name
            if assembly:
                raise ToolError(
                    f"this file is aligned to {assembly}, and lilr-wgs has no "
                    f"interval table for it.\n"
                    f"  {name} is {contigs[name]:,} bp; known assemblies are "
                    + ", ".join(f"{a} ({ln:,})"
                                for ln, a in sorted(CHR19_LENGTH.items())
                                if a in LOCI_BY_ASSEMBLY) + "\n"
                    "Chromosome 19 has the same name in every one of them, so "
                    "extracting anyway would return a different half-megabase "
                    "and call copy number on it. Lift the input over, or add a "
                    f"{assembly} interval table beside lilrwgs.loci_chm13."
                )
            raise ToolError(
                f"{name} is {contigs[name]:,} bp, which matches no assembly "
                "lilr-wgs knows.\n"
                f"  GRCh38: {58_617_616:,}   GRCh37: {59_128_983:,}\n"
                "If this is a patched or custom GRCh38, the LRC coordinates may "
                "still be right, but nothing here can confirm that."
            )
    raise ToolError(
        "no chromosome 19 in the header under any name lilr-wgs recognises.\n"
        f"  looked for: {', '.join(CHR19_ALIASES)}\n"
        f"  header has: {', '.join(sorted(contigs)[:6])} ...\n"
        "This is the whole input: LILR genes are on chr19q13.42 and nothing "
        "else in the file is used."
    )


def resolve_regions(contigs: dict[str, int], chr19: str, *,
                    include_alts: bool = True, loci_mod=loci,
                    ) -> tuple[list[str], bool]:
    """The slice intervals, spelled the way this header spells them.

    ``loci_mod`` is the table for the assembly the *input* is aligned to, which
    is not necessarily the one being realigned to: reads have to be found where
    they physically are before they can be moved. The coordinates are that
    table's and do not move; only the contig names are the input's.

    Alt contigs are included when present and dropped when not. CHM13 has none
    and needs none, so the flag is inert there.
    """
    alt_contigs = getattr(loci_mod, "ALT_CONTIGS", [])
    has_alts = any(c in contigs for c in alt_contigs)
    intervals = loci_mod.slice_intervals(include_alts=include_alts and has_alts)

    regions: list[str] = []
    for chrom, start, end in intervals:
        name = chr19 if chrom == loci_mod.CHROM else chrom
        if name in contigs:
            # Clip to the contig: a flanked interval can run past the end, and
            # samtools takes that as an error rather than as a truncation.
            regions.append(loci_mod.as_region(name, start,
                                              min(end, contigs[name])))
    return regions, has_alts


def assembly_of_reference(reference: str | os.PathLike) -> tuple[str, object]:
    """``(assembly, loci module)`` for a reference FASTA, from its ``.fai``.

    The target of a realignment is a FASTA, not a BAM header, so it cannot be
    identified the way :func:`source_assembly` identifies an input. The `.fai`
    carries the same fact — chromosome 19's length — and reading it costs
    nothing next to guessing from the filename, which is what a user renaming
    `chm13v2.0.fa` would break.
    """
    fai = Path(f"{reference}.fai")
    if not fai.exists():
        raise ToolError(
            f"{fai} not found; index the reference with `samtools faidx` so the "
            "assembly can be identified from it rather than from its name")
    lengths = {}
    for line in fai.read_text().splitlines():
        f = line.split("\t")
        if len(f) >= 2:
            lengths[f[0]] = int(f[1])
    assembly, _ = source_assembly(lengths)
    return assembly, LOCI_BY_ASSEMBLY[assembly]


def index_is_alt_aware(bwa_index: str | os.PathLike) -> bool:
    """Whether ``{index}.alt`` is beside the index, which is how bwa finds it.

    `bwa mem` takes no flag for this. It looks for the file, and silently aligns
    without ALT-awareness when it is absent — which produces a BAM where every
    MAPQ-20 window in the LRC reads near zero for every sample alike. Checked
    before a cohort rather than after it.
    """
    alt = Path(str(bwa_index) + ".alt")
    return alt.exists() and alt.stat().st_size > 0


def slice_to_fastq(
    sample: str,
    bam: str | os.PathLike,
    out_r1: str | os.PathLike,
    out_r2: str | os.PathLike,
    out_singletons: str | os.PathLike,
    *,
    threads: int = 4,
    tmpdir: str | os.PathLike | None = None,
    samtools: str = "samtools",
) -> tuple[int, int]:
    """The *whole* local slice back to FASTQ. Returns ``(n_pairs, n_singletons)``.

    Deliberately not :func:`lilrwgs.extract.to_fastq`, which restricts to the
    LILR cluster because it is feeding gene-level recruitment. Here the control
    loci are the point: they carry λ₁, and a realigned BAM without them has no
    baseline to divide a depth by.

    Singletons are written rather than discarded because they are realigned too.
    A pair whose mate fell outside the extracted intervals is still a read that
    contributed depth in the input alignment, and dropping it removes depth from
    the flanks of every interval — including the control loci, where the flank is
    1 kb and the effect is therefore largest relative to the interval.
    """
    require(samtools)
    bam = Path(bam).resolve()
    out_r1, out_r2 = Path(out_r1).resolve(), Path(out_r2).resolve()
    out_singletons = Path(out_singletons).resolve()
    for p in (out_r1, out_r2, out_singletons):
        p.parent.mkdir(parents=True, exist_ok=True)

    tmp = Path(tmpdir or os.environ.get("TMPDIR", "/tmp")) / f"lilrwgs_refq_{sample}"
    tmp.mkdir(parents=True, exist_ok=True)
    half = max(1, threads // 2)

    try:
        result = pipeline(
            [
                # -F SUPPLEMENTARY: a supplementary record is a fragment of a
                # read already in the file, and turning it into a FASTQ entry
                # would emit that read twice under one name. The alt-contig
                # evidence it represents is not lost — bwa regenerates it from
                # the full-length read at realignment, which is the point.
                [samtools, "view", "-u", "-F", str(SUPPLEMENTARY),
                 "-@", str(half), str(bam)],
                # collate, not `sort -n`: `samtools fastq` needs mates adjacent,
                # and grouping by name is all that requires — much cheaper than
                # a full sort.
                [samtools, "collate", "-u", "-O", "-@", str(half),
                 "-T", str(tmp / "collate"), "-"],
                # -n keeps read names as they came in, so a name in the realigned
                # BAM still matches the record it came from and a disagreement
                # can be traced back to one read.
                [samtools, "fastq", "-n", "-@", str(half),
                 "-1", str(out_r1), "-2", str(out_r2),
                 "-0", "/dev/null", "-s", str(out_singletons), "-"],
            ],
            cwd=str(tmp),
        )
        return _parse_fastq_counts(result.stderr)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def _parse_fastq_counts(stderr: str) -> tuple[int, int]:
    """Pull pair and singleton counts out of `samtools fastq` stderr.

    `processed` counts every read it saw, singletons included, so the pair count
    is ``(processed - singletons) / 2`` and not ``processed / 2``. Measured on
    HG00119: 167,082 processed and 2,362 singletons against 82,360 records in
    each of R1 and R2 — ``processed / 2`` gives 83,541, which is the pair count
    plus the singletons. The error is small, but it lands in the denominator of
    the singleton rate that decides whether to warn, so it is worth being right
    about rather than inheriting.

    The wording has shifted across samtools versions, so a miss returns zeros
    and the caller treats a zero-pair result as an error on its own terms rather
    than trusting this to have worked.
    """
    processed = singletons = 0
    for line in stderr.splitlines():
        low = line.lower()
        if "singleton" in low:
            for token in low.replace("[", " ").replace("]", " ").split():
                if token.isdigit():
                    singletons = int(token)
                    break
        elif "processed" in low and "read" in low:
            for token in low.split():
                if token.isdigit():
                    processed = int(token)
                    break
    return max(0, (processed - singletons) // 2), singletons


def align(
    sample: str,
    r1: str | os.PathLike,
    r2: str | os.PathLike,
    bwa_index: str | os.PathLike,
    out_bam: str | os.PathLike,
    *,
    singletons: str | os.PathLike | None = None,
    threads: int = 4,
    tmpdir: str | os.PathLike | None = None,
    bwa: str = "bwa",
    samtools: str = "samtools",
) -> int:
    """`bwa mem` the extracted reads back onto GRCh38. Returns record count.

    Two passes, merged. The paired one carries almost everything; the singleton
    one exists so that reads whose mate fell outside the extracted intervals
    still contribute the depth they contributed before, rather than being lost
    at the edge of every interval. bwa is given them separately because a
    single-end read in a paired run is not merely unpaired — it perturbs the
    insert-size estimate the paired reads are rescued with.
    """
    require(bwa, samtools)
    out_bam = Path(out_bam).resolve()
    out_bam.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(tmpdir or os.environ.get("TMPDIR", "/tmp")) / f"lilrwgs_bwa_{sample}"
    tmp.mkdir(parents=True, exist_ok=True)

    # Resolved here, not left to the caller. `bwa mem` runs with cwd set to the
    # scratch dir below, so a relative index base resolves against *that* and
    # bwa reports `fail locate index files` -- which names the index and says
    # nothing about the working directory, and sends you to check the index
    # files, which are fine. Same for the FASTQs.
    bwa_index = Path(bwa_index).resolve()
    r1, r2 = Path(r1).resolve(), Path(r2).resolve()
    if singletons:
        singletons = Path(singletons).resolve()

    # A read group, because the allele-calling path downstream needs one and
    # adding it here costs nothing. SM is what GATK keys a sample on.
    rg = f"@RG\\tID:{sample}\\tSM:{sample}\\tLB:{sample}\\tPL:ILLUMINA\\tPU:{sample}"

    try:
        parts: list[Path] = []
        jobs: list[tuple[str, list[str]]] = [
            ("paired", [str(r1), str(r2)]),
        ]
        if singletons and Path(singletons).exists() and Path(singletons).stat().st_size > 0:
            jobs.append(("single", [str(singletons)]))

        for label, inputs in jobs:
            part = tmp / f"{label}.bam"
            pipeline(
                [
                    [bwa, "mem", *BWA_ARGS, "-t", str(threads), "-R", rg,
                     str(bwa_index), *inputs],
                    [samtools, "sort", "-@", str(max(1, threads // 2)),
                     "-T", str(tmp / f"sort_{label}"), "-o", str(part), "-"],
                ],
                cwd=str(tmp),
            )
            parts.append(part)

        if len(parts) == 1:
            shutil.move(str(parts[0]), str(out_bam))
        else:
            run([samtools, "merge", "-f", "-@", str(threads),
                 str(out_bam), *[str(p) for p in parts]])
        run([samtools, "index", "-@", str(threads), str(out_bam)])
        return int(run([samtools, "view", "-c", str(out_bam)]).stdout.strip() or 0)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def realign_sample(
    sample: str,
    source: str,
    reference: str | os.PathLike,
    bwa_index: str | os.PathLike,
    out_bam: str | os.PathLike,
    *,
    threads: int = 4,
    tmpdir: str | os.PathLike | None = None,
    keep_slice: str | os.PathLike | None = None,
    index: str | os.PathLike | None = None,
    bwa: str = "bwa",
    samtools: str = "samtools",
) -> RealignStats:
    """One sample, from an arbitrary GRCh38 alignment to one this pipeline made.

    Args:
        source: CRAM or BAM, local path or ``https://`` URL.
        reference: the FASTA ``source`` was compressed against — needed to decode
            a CRAM, and unused for a BAM.
        bwa_index: the GRCh38 analysis-set index. ``{bwa_index}.alt`` must be
            beside it; without it the realignment is not ALT-aware and every
            MAPQ-20 window in the LRC reads near zero.
        out_bam: the realigned slice, in GRCh38 coordinates.
    """
    from .extract import slice_cram  # local: avoids a cycle at import time

    stats = RealignStats(sample=sample, source=str(source))
    stats.alt_aware = index_is_alt_aware(bwa_index)
    if not stats.alt_aware:
        stats.warnings.append(
            f"{bwa_index}.alt is missing, so the realignment is not ALT-aware "
            "and LILRA6/LILRB3 will be refused; run scripts/fetch_bwa_index.sh"
        )

    work = Path(tmpdir or os.environ.get("TMPDIR", "/tmp")) / f"lilrwgs_realign_{sample}"
    work.mkdir(parents=True, exist_ok=True)

    try:
        env = cram_env(reference)
        contigs = read_contigs(source, reference, samtools=samtools, env=env)
        stats.source_assembly, stats.chr19_name = source_assembly(contigs)
        regions, has_alts = resolve_regions(
            contigs, stats.chr19_name,
            loci_mod=LOCI_BY_ASSEMBLY[stats.source_assembly])
        stats.had_alt_contigs = has_alts

        started = time.time()
        slice_bam = work / f"{sample}.slice.bam"
        info = slice_cram(sample, source, reference, slice_bam, threads=threads,
                          tmpdir=str(work), regions=regions, index=index,
                          samtools=samtools)
        stats.n_source_records = info["n_records"]
        stats.slice_s = time.time() - started

        if keep_slice:
            keep = Path(keep_slice)
            keep.mkdir(parents=True, exist_ok=True)
            shutil.copy(slice_bam, keep / slice_bam.name)
            shutil.copy(str(slice_bam) + ".bai", keep / (slice_bam.name + ".bai"))

        r1, r2 = work / "lrc_R1.fq.gz", work / "lrc_R2.fq.gz"
        singles = work / "lrc_S.fq.gz"
        stats.n_pairs, stats.n_singletons = slice_to_fastq(
            sample, slice_bam, r1, r2, singles, threads=threads,
            tmpdir=str(work), samtools=samtools)
        if stats.n_pairs == 0:
            raise ToolError(
                f"{sample}: the slice held {stats.n_source_records} records but "
                "no read pairs survived collation — mates are not adjacent, or "
                "the records are all supplementary")
        if stats.n_singletons > 0.05 * stats.n_pairs:
            stats.warnings.append(
                f"singleton rate {stats.n_singletons / stats.n_pairs:.1%} is "
                "above 5%; mates are falling outside the extracted intervals "
                "more often than the insert size explains")

        started = time.time()
        stats.n_realigned = align(
            sample, r1, r2, bwa_index, out_bam, singletons=singles,
            threads=threads, tmpdir=str(work), bwa=bwa, samtools=samtools)
        stats.align_s = time.time() - started

        if stats.n_realigned == 0:
            raise ToolError(f"{sample}: bwa produced an empty BAM from "
                            f"{stats.n_pairs} pairs")
        return stats
    finally:
        shutil.rmtree(work, ignore_errors=True)
