"""CRAM -> a local slice -> LRC read pairs as FASTQ.

The front end. Where `lilr-genotyper` started from a whole readset and ran 11
full-readset alignments to find LILR reads, srWGS input is already aligned, so
the same job is one indexed slice: the reads are where the reference says they
are.

That is only true if the slice covers everywhere a LILR read could have landed,
which is why :func:`lilrwgs.loci.extraction_regions` returns the primary cluster
*and* the alt-contig intervals carrying LILRA3. A read from a LILRA3-bearing
chromosome has its only primary alignment there, so a slice of chr19 alone would
drop LILRA3 entirely and do it silently — every sample would look like a
deletion homozygote.

Two steps, not one. :func:`slice_cram` makes a single remote pass and writes a
local, indexed BAM covering the union of everything any later stage needs;
:func:`to_fastq` and the coverage model then read that. The split exists because
the coverage model's control loci sit deliberately *outside* the LILR cluster,
so a slice shaped for extraction would not contain them, and a CRAM read over
HTTPS is the expensive part of a sample — at 2,504 samples, doing it twice is
not a rounding error. The local BAM is also a resumption point: everything after
it is cheap to redo.

Nothing is staged beyond that slice. htslib reads the CRAM over HTTPS using the
`.crai`, fetching only containers that overlap the requested regions, so a sample
costs tens of megabytes against a CRAM of tens of gigabytes. This needs a
samtools built with libcurl; :func:`check_remote_support` says so before a cohort
discovers it the hard way.
"""

from __future__ import annotations

import os
import shutil
import time
from dataclasses import dataclass, field
from pathlib import Path

from . import loci
from .shell import ToolError, Result, pipeline, require, run

# Records excluded at the slice: secondary (0x100), QC-fail (0x200),
# duplicate (0x400), supplementary (0x800).
#
# Duplicate handling is one of the sharpest differences from the capture
# pipeline, which skips MarkDuplicates outright: in targeted capture, read pairs
# legitimately share start coordinates because probe geometry drives where reads
# begin, and marking flagged 70-80% of pairs and crushed depth to nothing. In
# WGS a shared start coordinate means what it usually means. The 1000 Genomes
# CRAMs arrive duplicate-marked by NYGC, so the flags are honoured rather than
# recomputed — and recomputing them on a 550 kb slice would be wrong anyway,
# since duplicate detection needs the whole library to judge.
EXCLUDE_FLAGS = 0xF00


@dataclass
class ExtractionStats:
    """What came out, and enough context to know whether to believe it."""

    sample: str
    source: str
    regions: list[str]
    n_pairs: int = 0
    n_singletons: int = 0
    elapsed_s: float = 0.0
    r1: str = ""
    r2: str = ""
    singletons: str = ""
    warnings: list[str] = field(default_factory=list)

    def as_row(self) -> dict:
        return {
            "sample": self.sample,
            "n_pairs": self.n_pairs,
            "n_singletons": self.n_singletons,
            "singleton_rate": round(self.n_singletons / max(self.n_pairs, 1), 5),
            "elapsed_s": round(self.elapsed_s, 1),
            "n_regions": len(self.regions),
            "warnings": ";".join(self.warnings),
        }


def cram_env(reference: str | os.PathLike | None = None,
             base: dict | None = None) -> dict:
    """Environment for an htslib process reading a CRAM.

    ``REF_PATH`` is pinned rather than left to htslib's default, which points at
    the EBI CRAM reference registry. On a host that cannot reach it — or one that
    can, but slowly, 2,504 times — a silent outbound fetch per container hangs
    instead of failing, and the symptom is a job that looks merely slow. Setting
    it to a path that does not resolve forces htslib to use the ``-T`` reference
    we passed it, and to fail loudly if that reference is wrong.
    """
    env = dict(base if base is not None else os.environ)
    env["REF_PATH"] = str(reference) if reference else "/dev/null/no_ref_registry"
    return env


def check_remote_support(samtools: str = "samtools") -> None:
    """Fail if this samtools cannot read a URL.

    Without libcurl, `samtools view https://...` fails with a bare "fail to open
    file", which reads exactly like a typo in a path and sends people looking in
    the wrong place. The distribution matters: the environment-module samtools on
    at least one HPC site is built ``--without-libcurl``, so this is a real
    configuration, not a hypothetical.
    """
    require(samtools)
    version = run([samtools, "--version"]).stdout
    # `samtools --version` prints a Features line for samtools itself and
    # another for htslib. libcurl belongs to htslib, and on a build where the
    # two differ the samtools line reads `build=configure curses=yes` with no
    # mention of it — so checking only the first line reports a perfectly good
    # installation as broken.
    features = [ln.strip() for ln in version.splitlines()
                if ln.strip().startswith("Features:")]
    if not features:
        raise ToolError(
            f"could not parse `{samtools} --version` for a Features line")
    if not any("libcurl=yes" in line for line in features):
        raise ToolError(
            f"{shutil.which(samtools)} is built without libcurl.\n"
            + "\n".join(f"  {line}" for line in features) + "\n"
            "It cannot read CRAMs over HTTPS, which is how this pipeline reads "
            "1000 Genomes data. Use the conda environment (environment.yml) "
            "rather than a system or module build, or stage CRAMs locally and "
            "point the manifest at the files."
        )


def slice_cram(
    sample: str,
    cram: str,
    reference: str | os.PathLike,
    out_bam: str | os.PathLike,
    *,
    threads: int = 4,
    tmpdir: str | os.PathLike | None = None,
    include_alts: bool = True,
    samtools: str = "samtools",
) -> dict:
    """The one remote pass. Fetch every region any stage needs into a local BAM.

    Args:
        cram: local path or ``https://`` URL.
        reference: the FASTA the CRAM was compressed against. For 1000 Genomes
            that is ``GRCh38_full_analysis_set_plus_decoy_hla.fa``. A *different*
            GRCh38 will decode most bases correctly and corrupt the rest, which
            is considerably worse than failing.
        include_alts: include the alt-contig intervals carrying LILRA3. Off only
            makes sense for a CRAM aligned to a primary-only reference, where
            LILRA3 is then not measurable by depth at all.

    Returns:
        A dict of counts and timings, including whether each contig class was
        present in the CRAM header.
    """
    require(samtools)
    out_bam = Path(out_bam).resolve()
    out_bam.parent.mkdir(parents=True, exist_ok=True)
    # The stages below run in a scratch cwd, so every path handed to a tool has
    # to be absolute. A relative reference resolves against the scratch dir and
    # fails with "Failed to open reference file", which points at the reference
    # rather than at the cwd and is a genuinely confusing place to start looking.
    reference = Path(reference).resolve()
    tmp = Path(tmpdir or os.environ.get("TMPDIR", "/tmp")) / f"lilrwgs_slice_{sample}"
    tmp.mkdir(parents=True, exist_ok=True)
    env = cram_env(reference)
    started = time.time()

    # Which contigs the CRAM actually has. Asking up front separates "this CRAM
    # was aligned to a primary-only reference" — a real, reportable fact about
    # the data — from "samtools failed", which a bare error on an unknown region
    # would conflate. Both produce no LILRA3; only one is a bug.
    header = run([samtools, "view", "-H", "-T", str(reference), cram], env=env).stdout
    contigs = {part[3:] for ln in header.splitlines() if ln.startswith("@SQ")
               for part in ln.split("\t") if part.startswith("SN:")}

    has_alts = any(c in contigs for c in loci.ALT_CONTIGS)
    if include_alts and not has_alts:
        include_alts = False

    regions = loci.slice_regions(include_alts=include_alts)
    regions = [r for r in regions if r.split(":")[0] in contigs]
    if not regions:
        raise ToolError(
            f"{sample}: none of the requested contigs are in the CRAM header.\n"
            f"  looked for: {loci.CHROM}, {loci.ALT_CONTIGS[0]}, ...\n"
            f"  header has: {', '.join(sorted(contigs)[:6])} ...\n"
            "This is usually Ensembl-style naming ('19' not 'chr19'); "
            "lilr-wgs expects the UCSC-style names the 1000 Genomes CRAMs use."
        )

    try:
        pipeline(
            [
                # -M: with several regions, samtools emits a read once per
                # region it overlaps unless told otherwise. The regions are
                # merged upstream, but overlap here doubles depth rather than
                # erroring, so both guards are kept.
                [samtools, "view", "-u", "-M", "-T", str(reference),
                 "-F", str(EXCLUDE_FLAGS),
                 "-@", str(max(1, threads // 2)), cram, *regions],
                [samtools, "sort", "-@", str(max(1, threads // 2)),
                 "-T", str(tmp / "sort"), "-o", str(out_bam), "-"],
            ],
            env=env, cwd=str(tmp),
        )
        run([samtools, "index", "-@", str(threads), str(out_bam)])
        n_records = int(run([samtools, "view", "-c", str(out_bam)]).stdout.strip() or 0)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    if n_records == 0:
        raise ToolError(
            f"{sample}: the slice is empty.\n"
            f"  source: {cram}\n"
            f"  regions: {len(regions)}\n"
            "The contig names matched, so this is most likely a reference "
            "mismatch on -T, or a CRAM whose index does not cover these regions."
        )

    return {
        "sample": sample,
        "source": cram,
        "bam": str(out_bam),
        "n_records": n_records,
        "n_regions": len(regions),
        "has_alt_contigs": has_alts,
        "elapsed_s": round(time.time() - started, 1),
    }


def to_fastq(
    sample: str,
    bam: str | os.PathLike,
    out_r1: str | os.PathLike,
    out_r2: str | os.PathLike,
    *,
    out_singletons: str | os.PathLike | None = None,
    threads: int = 4,
    tmpdir: str | os.PathLike | None = None,
    include_alts: bool = True,
    samtools: str = "samtools",
) -> ExtractionStats:
    """Turn the LILR-bearing part of a local slice into paired FASTQ.

    Restricted to :func:`lilrwgs.loci.extraction_regions` — the slice also holds
    the coverage model's control loci, and those reads have no business in a
    gene-assignment FASTQ.

    Returns:
        ExtractionStats including a singleton count. Singletons are pairs whose
        mate fell outside the extracted regions; with ~65 kb of flank and a
        450 bp insert they should be a rounding error, and a high rate means
        either the flanks are wrong or the library is not what we think it is —
        so the rate is reported rather than quietly discarded.
    """
    require(samtools)
    # Absolute, for the same reason as in slice_cram: these stages run in a
    # scratch cwd.
    bam = Path(bam).resolve()
    out_r1, out_r2 = Path(out_r1).resolve(), Path(out_r2).resolve()
    out_r1.parent.mkdir(parents=True, exist_ok=True)
    out_r2.parent.mkdir(parents=True, exist_ok=True)

    regions = loci.extraction_regions(include_alts=include_alts)
    tmp = Path(tmpdir or os.environ.get("TMPDIR", "/tmp")) / f"lilrwgs_fastq_{sample}"
    tmp.mkdir(parents=True, exist_ok=True)
    singles = Path(out_singletons) if out_singletons else tmp / "singletons.fq.gz"
    half = max(1, threads // 2)
    started = time.time()

    stages = [
        # -u: uncompressed BAM between stages; the pipe is local and CPU time is
        # worth more here than the bytes.
        [samtools, "view", "-u", "-M", "-@", str(half), str(bam), *regions],
        # The slice is coordinate-sorted and `samtools fastq` needs mates
        # adjacent. collate, not `sort -n`: it groups by name without a full
        # sort, which is all that is needed and much cheaper.
        [samtools, "collate", "-u", "-O", "-@", str(half), "-T", str(tmp / "collate"), "-"],
        # -n keeps read names exactly as they came in, so a name in the FASTQ
        # still matches the CRAM record it came from and a disagreement
        # downstream can be traced back to one read.
        [samtools, "fastq", "-n", "-@", str(half),
         "-1", str(out_r1), "-2", str(out_r2),
         "-0", "/dev/null", "-s", str(singles), "-"],
    ]

    result: Result = pipeline(stages, cwd=str(tmp))

    stats = ExtractionStats(
        sample=sample, source=str(bam), regions=regions,
        elapsed_s=time.time() - started,
        r1=str(out_r1), r2=str(out_r2), singletons=str(singles),
    )
    stats.n_pairs, stats.n_singletons = _parse_fastq_counts(result.stderr)

    if stats.n_pairs == 0:
        raise ToolError(
            f"{sample}: no read pairs in the LILR regions of {bam}.\n"
            "The slice was non-empty, so the reads are in the control loci and "
            "not the cluster — check that the extraction regions and the slice "
            "were built for the same assembly."
        )
    if stats.n_singletons > 0.05 * stats.n_pairs:
        stats.warnings.append(
            f"singleton rate {stats.n_singletons / stats.n_pairs:.1%} is above 5%; "
            "mates are falling outside the extracted regions more often than the "
            "insert size explains"
        )

    shutil.rmtree(tmp, ignore_errors=True)
    return stats


def _parse_fastq_counts(stderr: str) -> tuple[int, int]:
    """Pull pair and singleton counts out of `samtools fastq` stderr.

    It reports e.g. "processed 41234 reads" and "discarded N singletons". The
    format has shifted across samtools versions, so a miss returns zeros and the
    caller treats a zero-pair result as an error on its own terms rather than
    trusting this to have worked.
    """
    pairs = singletons = 0
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
                    pairs = int(token) // 2
                    break
    return pairs, singletons


def main(argv: list[str] | None = None) -> int:
    import argparse
    import json

    p = argparse.ArgumentParser(
        description="Slice the LRC out of a CRAM and write paired FASTQ")
    p.add_argument("sample")
    p.add_argument("cram", help="local path or https:// URL")
    p.add_argument("reference", help="FASTA the CRAM was compressed against")
    p.add_argument("out_bam", help="the local slice, kept as a resumption point")
    p.add_argument("out_r1")
    p.add_argument("out_r2")
    p.add_argument("--singletons")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--tmpdir")
    p.add_argument("--no-alts", action="store_true",
                   help="skip the LILRA3 alt-contig intervals, making LILRA3 "
                        "unmeasurable by depth")
    p.add_argument("--stats", help="write the stats row here as JSON")
    args = p.parse_args(argv)

    if args.cram.startswith(("http://", "https://", "s3://", "gs://")):
        check_remote_support()

    slice_info = slice_cram(
        args.sample, args.cram, args.reference, args.out_bam,
        threads=args.threads, tmpdir=args.tmpdir, include_alts=not args.no_alts,
    )
    if not slice_info["has_alt_contigs"] and not args.no_alts:
        print(f"warning: {args.sample}: the CRAM's reference has no LRC alt "
              "contigs, so LILRA3 has no depth window; only the junction assay "
              "will be available")

    stats = to_fastq(
        args.sample, args.out_bam, args.out_r1, args.out_r2,
        out_singletons=args.singletons, threads=args.threads, tmpdir=args.tmpdir,
        include_alts=slice_info["has_alt_contigs"] and not args.no_alts,
    )
    row = {**slice_info, **stats.as_row()}
    if args.stats:
        Path(args.stats).parent.mkdir(parents=True, exist_ok=True)
        Path(args.stats).write_text(json.dumps(row, indent=2))
    for warning in stats.warnings:
        print(f"warning: {args.sample}: {warning}")
    print(f"{args.sample}: {row['n_records']} records sliced in "
          f"{slice_info['elapsed_s']}s, {row['n_pairs']} pairs written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
