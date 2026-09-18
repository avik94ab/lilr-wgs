"""Measuring what one haploid copy looks like in this sample.

Everything in :mod:`lilrwgs.depth_model` is expressed in units of λ₁, the depth a
single haploid copy of unique sequence yields. This module measures it, along
with the spread around it and two things that decide whether the measurement can
be trusted at all.

λ₁ comes from control loci rather than from the LILR genes themselves, and that
is the structural advantage srWGS has over targeted capture. In a capture panel
the only available reference for "what does one copy look like?" is other genes
in the same panel — so copy number is inherently relative, a cohort has to be
called as one batch, and a panel-wide capture efficiency shift is invisible.
Here the controls are ordinary genomic sequence a few hundred kilobases away,
measured in the same CRAM, and a single sample can be called on its own.

Two diagnostics decide whether the MAPQ-20 numbers mean anything:

**Was the alignment ALT-aware?** GRCh38 places nine LRC alt haplotypes at one
primary interval that contains every window used here. Given the ``.alt`` file,
bwa computes MAPQ over primary hits alone and a read in an LRC window keeps it;
without it, the same read has nine more equally good placements and comes back at
MAPQ 0. In that case every MAPQ-20 count in the cluster reads near zero — for
everyone, targets and controls alike — which is indistinguishable from a cohort
of homozygous deletions unless something outside the placement is measured too.
That is what :data:`lilrwgs.loci.OUTSIDE_LOCI` is for.

**Is the region being sampled normally at all?** A sample whose control depth is
wildly off the cohort, or whose GC curve is steep, is not one whose LILR calls
should carry the same confidence as everyone else's.
"""

from __future__ import annotations

import math
import os
import statistics
from dataclasses import dataclass, field
from pathlib import Path

from . import loci
from .depth_model import estimate_dispersion
from .shell import ToolError, require, run

# Bases below this quality are not counted. 13 is samtools' own default for
# mpileup and the same threshold the depth is modelled on, so the number that
# reaches the model is the number a caller would see.
MIN_BASEQ = 13

# GC is computed in windows this wide around each position. 100 bp is the read
# length, which is the scale at which GC actually affects whether a fragment was
# amplified and sequenced.
GC_WINDOW = 100

# A GC bin needs this many bases behind it before its correction is used. Below
# that the bin's median is noise, and a noisy multiplier on λ₁ is worse than no
# correction — it would move thresholds for reasons unrelated to the sample.
MIN_BASES_PER_GC_BIN = 2_000

# Verdict thresholds, from the LRC-vs-outside comparison.
MIN_Q20_LRC = 0.30          # below this, MAPQ-20 windows inside the LRC are dead
MAX_DILUTION = 1.50         # above this, reads are scattering across alt haplotypes


# Genes used to calibrate the recruitment path. Separable from their neighbours
# (so no shared-block inflation) and copy-number stable at 2 in >98% of
# pangenome samples. LILRB1 was an anchor in the predecessor's CN caller but is
# excluded here: it shares a block with LILRB4, so its depth carries both genes
# and would bias the factor upward.
EFFICIENCY_ANCHORS = ["LILRA1", "LILRA2", "LILRB2", "LILRB5"]


def recruitment_efficiency(anchor_depths: dict[str, float], lambda1: float,
                           ) -> tuple[float, int]:
    """How much depth survives the path from CRAM to per-gene alignment.

    lambda_1 is measured on the CRAM slice, but the depth a threshold is applied
    to is measured much later — after extraction to FASTQ, recruitment against a
    pangenome panel, cross-map arbitration, and realignment to a single per-locus
    reference. Every one of those steps drops reads, so the two numbers are in
    different units and comparing them directly makes every gene look
    under-covered.

    Measured on HG00096: the CRAM says lambda_1 = 19.1, so a diploid locus should
    carry 38.2, and LILRB1's realigned median was 31 — about 19% down. Applied
    without this correction, that shortfall becomes 27% of the gene falling below
    a floor it should have cleared, which is the predecessor's failure mode
    reappearing by a different route.

    This does not reintroduce a relative baseline. The *copy number* still comes
    from the CRAM-level measurement against external controls; what is estimated
    here is the efficiency of a fixed pipeline path, which is a property of the
    software and not of the sample's biology. The anchors are separable,
    copy-stable genes, so their expected depth is 2 * lambda_1 by construction.

    Returns:
        (factor, n_anchors). The factor is clamped to a sane band: a value far
        from 1 means the recruitment path is broken rather than merely lossy, and
        silently scaling by it would hide that.
    """
    usable = [d for d in anchor_depths.values() if d > 0]
    if not usable or lambda1 <= 0:
        return 1.0, 0
    observed = statistics.median(usable)
    factor = observed / (2.0 * lambda1)
    return min(1.5, max(0.3, factor)), len(usable)


@dataclass
class ControlMeasurement:
    """Per-base depth over one control locus, at both MAPQ floors."""

    name: str
    chrom: str
    start: int
    end: int
    inside_placement: bool
    mean_q0: float = 0.0
    mean_q20: float = 0.0
    var_q20: float = 0.0
    n_bases: int = 0

    @property
    def q20_retention(self) -> float:
        """Fraction of depth surviving the MAPQ floor. ~1 where a floor works."""
        return self.mean_q20 / self.mean_q0 if self.mean_q0 > 0 else 0.0


@dataclass
class CoverageModel:
    """What one haploid copy looks like here, and whether to believe it."""

    sample: str
    lambda1: float = 0.0
    dispersion: float = math.inf
    lambda1_outside: float = 0.0
    q20_lrc: float = 0.0
    q20_outside: float = 0.0
    dilution: float = 0.0
    alt_verdict: str = "unknown"
    # Depth surviving the CRAM -> FASTQ -> panel -> locus-reference path.
    # 1.0 until measured; see recruitment_efficiency().
    efficiency: float = 1.0
    n_efficiency_anchors: int = 0
    gc_correction: dict[int, float] = field(default_factory=dict)
    controls: list[ControlMeasurement] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def usable_mapq20(self) -> bool:
        """Whether MAPQ-20 measurements inside the LRC mean anything.

        When they do not, the honest report is "not measured". Dividing a dead
        MAPQ-20 count by a live baseline produces a clean, confident zero, and a
        cohort of zeroes at LILRA6 looks exactly like a cohort of deletion
        homozygotes — a wrong answer that arrives looking like a right one.
        """
        return self.q20_lrc >= MIN_Q20_LRC

    def lambda_at(self, gc: float | None = None) -> float:
        """λ₁ at a position of the given GC fraction, in per-gene-BAM units.

        The efficiency factor is applied here rather than left to callers,
        because every consumer of λ₁ works on a realigned per-gene BAM and would
        otherwise have to remember to apply it — and forgetting would look like
        a gene with poor coverage rather than like a missing correction.
        """
        lam = self.lambda1 * self.efficiency
        if gc is None or not self.gc_correction:
            return lam
        return lam * self.gc_correction.get(gc_bin(gc), 1.0)

    def as_row(self) -> dict:
        return {
            "sample": self.sample,
            "lambda1": round(self.lambda1, 3),
            "dispersion": (round(self.dispersion, 2)
                           if math.isfinite(self.dispersion) else "inf"),
            "lambda1_outside": round(self.lambda1_outside, 3),
            "q20_lrc": round(self.q20_lrc, 4),
            "q20_outside": round(self.q20_outside, 4),
            "dilution": round(self.dilution, 4),
            "alt_verdict": self.alt_verdict,
            "efficiency": round(self.efficiency, 4),
            "n_efficiency_anchors": self.n_efficiency_anchors,
            "usable_mapq20": self.usable_mapq20,
            "n_controls": len(self.controls),
            "warnings": ";".join(self.warnings),
        }


def gc_bin(gc: float) -> int:
    """GC fraction -> bin index in 5% steps."""
    return min(20, max(0, int(round(gc * 100 / 5))))


def alignment_verdict(q20_lrc: float, q20_outside: float, dilution: float) -> str:
    """Name what the diagnostics are saying, in words.

    Three numbers, and their combinations mean different things:

    - both retentions healthy, dilution ~1: an ALT-aware alignment, everything
      works.
    - outside healthy, LRC collapsed: the alt contigs are eating MAPQ inside the
      placement. LILRA3 survives on the outside baseline; LILRA6 and LILRB3 do
      not, and must be reported as not measured rather than as zero.
    - both collapsed: MAPQ is degraded genome-wide, which is a property of the
      sample or the aligner and not of the alt contigs.
    """
    if q20_outside < MIN_Q20_LRC:
        return "mapq_degraded_everywhere"
    if q20_lrc < MIN_Q20_LRC:
        return "not_alt_aware"
    if dilution > MAX_DILUTION:
        return "alt_diluted"
    return "alt_aware"


def fit_lambda(controls: list[ControlMeasurement]) -> tuple[float, float]:
    """λ₁ and the per-haploid dispersion from a set of controls.

    The controls sit at copy number 2, so λ₁ is half their depth. Dispersion is
    additive across copies — NB(μ, r) + NB(μ, r) is NB(2μ, 2r) — so the
    per-haploid ``r`` is half the value fitted at the controls' own depth.

    The median is taken across loci rather than the mean: one control overlapping
    an unannotated CNV in one sample should move the estimate by nothing, and
    with four loci the mean would move by a quarter of the error.
    """
    usable = [c for c in controls if c.n_bases > 0 and c.mean_q20 > 0]
    if not usable:
        return 0.0, math.inf

    mean_q20 = statistics.median(c.mean_q20 for c in usable)
    lambda1 = mean_q20 / 2.0

    # Dispersion from the control with the most bases behind it, rather than a
    # median of per-locus estimates: r is a ratio of moments and medians of
    # ratios do not compose the way medians of depths do.
    anchor = max(usable, key=lambda c: c.n_bases)
    r_diploid = estimate_dispersion(anchor.mean_q20, anchor.var_q20, anchor.n_bases)
    dispersion = r_diploid / 2.0 if math.isfinite(r_diploid) else math.inf
    return lambda1, dispersion


def fit_gc_correction(gc_by_base: list[float], depth_by_base: list[int],
                      ) -> dict[int, float]:
    """Multiplicative depth correction by GC bin.

    A 30x PCR-free library still varies by 10-30% across the GC range, and the
    LILR genes are not at the genomic mean GC — so an uncorrected λ₁ is
    systematically wrong at the loci this pipeline cares about, in a direction
    that depends on the gene. Bins with too little support are omitted and
    :meth:`CoverageModel.lambda_at` falls back to 1.0 for them.

    Two different statistics, deliberately:

    - each bin's own level is a **median**, so one pile-up inside a bin does not
      drag that bin's correction;
    - the normaliser is the **mean** over all control bases, because that is what
      :func:`fit_lambda` averages to get λ₁. Normalising by the median instead
      would make ``lambda_at`` disagree with ``lambda1`` by whatever the
      mean-median gap happens to be — a silent, GC-dependent bias in every
      threshold, which is precisely the kind of error this module exists to
      remove.
    """
    if not depth_by_base:
        return {}

    buckets: dict[int, list[int]] = {}
    for gc, depth in zip(gc_by_base, depth_by_base):
        buckets.setdefault(gc_bin(gc), []).append(depth)

    overall = statistics.fmean(depth_by_base)
    if overall <= 0:
        return {}

    return {
        b: statistics.median(depths) / overall
        for b, depths in buckets.items()
        if len(depths) >= MIN_BASES_PER_GC_BIN and statistics.median(depths) > 0
    }


def _depth_by_region(bam: str | os.PathLike, regions: list[tuple[str, int, int]],
                     mapq: int, *, reference: str | os.PathLike | None = None,
                     samtools: str = "samtools",
                     ) -> dict[str, list[int]]:
    """Per-base depth for each region, keyed by ``chrom:start-end``.

    ``-a`` matters: without it samtools omits zero-depth positions, and a region
    that is half uncovered would report the mean of its covered half. That is the
    difference between "this control looks normal" and "half of this control is
    missing", which is exactly the kind of thing the model must not be blind to.
    """
    require(samtools)
    bed_lines = "".join(f"{c}\t{s}\t{e}\n" for c, s, e in regions)
    # -G 0x800 excludes supplementary. The slice keeps them because LILRA3's
    # evidence lives there, but lambda_1 must be measured the way it always was:
    # at unique control sequence a supplementary record is a misplaced fragment,
    # not extra coverage, and counting it would inflate the baseline that every
    # threshold is expressed in.
    cmd = [samtools, "depth", "-a", "-q", str(MIN_BASEQ), "-Q", str(mapq),
           "-G", "0x800", "-b", "/dev/stdin"]
    if reference:
        cmd += ["--reference", str(reference)]
    cmd.append(str(bam))

    out = run(cmd, text_input=bed_lines).stdout

    depths: dict[str, list[int]] = {}
    lookup = {(c, s, e): f"{c}:{s}-{e}" for c, s, e in regions}
    by_chrom: dict[str, list[tuple[int, int, str]]] = {}
    for (c, s, e), key in lookup.items():
        by_chrom.setdefault(c, []).append((s, e, key))
        depths[key] = [0] * (e - s)

    for line in out.splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        chrom, pos, depth = parts[0], int(parts[1]) - 1, int(parts[2])
        for s, e, key in by_chrom.get(chrom, ()):
            if s <= pos < e:
                depths[key][pos - s] = depth
                break
    return depths


def measure(
    sample: str,
    bam: str | os.PathLike,
    *,
    reference: str | os.PathLike | None = None,
    samtools: str = "samtools",
    with_gc: bool = True,
    loci_mod=loci,
) -> CoverageModel:
    """Build the coverage model for one sample from its local slice.

    Args:
        bam: the local slice from :func:`lilrwgs.extract.slice_cram`. It must
            contain the control loci, which is why the slice is built from
            :func:`lilrwgs.loci.slice_regions` rather than the extraction
            regions alone.
        reference: needed only to compute the GC curve.
        loci_mod: the coordinate table `bam` is expressed in —
            :mod:`lilrwgs.loci` for GRCh38, :mod:`lilrwgs.loci_chm13` for
            CHM13v2.0. Passing the wrong one does not fail: both assemblies call
            the chromosome chr19 and both coordinates exist, so the controls are
            simply read ~3 Mb from where they live and λ₁ comes back as
            whatever sequence happens to be there. The caller that opens the BAM
            knows which reference made it; this argument is how it says so.
    """
    model = CoverageModel(sample=sample)

    regions = [(c.chrom, c.start, c.end) for c in loci_mod.ALL_CONTROLS]
    q0 = _depth_by_region(bam, regions, loci_mod.MAPQ_ANY,
                          reference=reference, samtools=samtools)
    q20 = _depth_by_region(bam, regions, loci_mod.MAPQ_STRICT,
                           reference=reference, samtools=samtools)

    for control in loci_mod.ALL_CONTROLS:
        key = f"{control.chrom}:{control.start}-{control.end}"
        d0, d20 = q0.get(key, []), q20.get(key, [])
        if not d20:
            continue
        model.controls.append(ControlMeasurement(
            name=control.name, chrom=control.chrom,
            start=control.start, end=control.end,
            inside_placement=control.inside_placement,
            mean_q0=statistics.fmean(d0) if d0 else 0.0,
            mean_q20=statistics.fmean(d20),
            var_q20=statistics.variance(d20) if len(d20) > 1 else 0.0,
            n_bases=len(d20),
        ))

    if not model.controls:
        raise ToolError(
            f"{sample}: no control loci could be measured in {bam}.\n"
            "The slice is missing the control regions — it was probably built "
            "from extraction_regions() rather than slice_regions()."
        )

    inside = [c for c in model.controls if c.inside_placement]
    outside = [c for c in model.controls if not c.inside_placement]

    model.lambda1, model.dispersion = fit_lambda(inside or model.controls)
    if outside:
        model.lambda1_outside, _ = fit_lambda(outside)

    model.q20_lrc = statistics.median([c.q20_retention for c in inside]) if inside else 0.0
    model.q20_outside = (statistics.median([c.q20_retention for c in outside])
                         if outside else 0.0)

    # Dilution: MAPQ-0 depth per base outside the placement over inside it. ~1
    # in an ALT-aware alignment; well above 1 when reads that belong in the LRC
    # are being scattered across the alt haplotypes, which rescales every MAPQ-0
    # ratio normalised on an LRC control.
    mean_in = statistics.median([c.mean_q0 for c in inside]) if inside else 0.0
    mean_out = statistics.median([c.mean_q0 for c in outside]) if outside else 0.0
    model.dilution = mean_out / mean_in if mean_in > 0 else 0.0
    model.alt_verdict = alignment_verdict(model.q20_lrc, model.q20_outside,
                                          model.dilution)

    if with_gc and reference:
        model.gc_correction = _gc_from_controls(q20, reference)

    _add_warnings(model)
    return model


def _gc_from_controls(depths_by_region: dict[str, list[int]],
                      reference: str | os.PathLike) -> dict[int, float]:
    """GC curve from the control loci's per-base depth."""
    try:
        import pysam
    except ImportError:
        return {}

    gc_vals: list[float] = []
    depth_vals: list[int] = []
    half = GC_WINDOW // 2
    with pysam.FastaFile(str(reference)) as fa:
        for key, depths in depths_by_region.items():
            chrom, span = key.rsplit(":", 1)
            start = int(span.split("-")[0])
            seq = fa.fetch(chrom, max(0, start - half),
                           start + len(depths) + half).upper()
            for i, depth in enumerate(depths):
                window = seq[i:i + GC_WINDOW]
                if len(window) < GC_WINDOW:
                    continue
                acgt = window.count("A") + window.count("C") + \
                    window.count("G") + window.count("T")
                if acgt < GC_WINDOW * 0.9:      # N-rich, GC is not meaningful
                    continue
                gc_vals.append((window.count("G") + window.count("C")) / acgt)
                depth_vals.append(depth)
    return fit_gc_correction(gc_vals, depth_vals)


def _add_warnings(model: CoverageModel) -> None:
    if model.lambda1 <= 0:
        model.warnings.append("lambda1 is zero: no usable control depth")
    elif model.lambda1 < 8:
        model.warnings.append(
            f"lambda1 {model.lambda1:.1f} is low for a nominally 30x library; "
            "expect a high uncallable fraction"
        )
    if model.alt_verdict == "not_alt_aware":
        model.warnings.append(
            "MAPQ 20 is dead inside the LRC but alive outside it: the alignment "
            "was not ALT-aware. LILRA6 and LILRB3 are not measurable by depth; "
            "they will be reported as not measured rather than as zero"
        )
    elif model.alt_verdict == "mapq_degraded_everywhere":
        model.warnings.append(
            "MAPQ 20 retention is low outside the alt placement too, so this is "
            "not an alt-contig problem — treat every MAPQ-filtered number here "
            "with suspicion"
        )
    elif model.alt_verdict == "alt_diluted":
        model.warnings.append(
            f"dilution {model.dilution:.2f}: MAPQ-0 depth is markedly lower "
            "inside the alt placement than outside, so reads are scattering "
            "across the alt haplotypes"
        )


def main(argv: list[str] | None = None) -> int:
    import argparse
    import json

    p = argparse.ArgumentParser(description="Fit the per-sample coverage model")
    p.add_argument("sample")
    p.add_argument("bam", help="local slice from `lilrwgs.extract`")
    p.add_argument("-o", "--output", help="write the model as JSON")
    p.add_argument("--reference", help="needed for the GC correction")
    p.add_argument("--no-gc", action="store_true")
    args = p.parse_args(argv)

    model = measure(args.sample, args.bam, reference=args.reference,
                    with_gc=not args.no_gc)
    row = model.as_row()
    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        Path(args.output).write_text(json.dumps({
            **row,
            "gc_correction": model.gc_correction,
            "controls": [vars(c) for c in model.controls],
        }, indent=2))
    for w in model.warnings:
        print(f"warning: {args.sample}: {w}")
    print(f"{args.sample}: lambda1={row['lambda1']} dispersion={row['dispersion']} "
          f"verdict={row['alt_verdict']} (q20_lrc={row['q20_lrc']}, "
          f"q20_out={row['q20_outside']}, dilution={row['dilution']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
