"""Copy number for LILRA3, LILRA6 and LILRB3.

The three variable genes, each difficult in its own way, each measured by more
than one route so that the routes can be made to disagree.

**LILRA6 and LILRB3** are read by depth over their paralogue-unique windows —
2,900 bp and 1,881 bp at the 3' ends, where an activating receptor and an
inhibitory one genuinely differ — divided by λ₁. Those windows only exist at
MAPQ 20, so the whole measurement is conditional on the alignment having been
ALT-aware, and :attr:`CoverageModel.usable_mapq20` gates it.

**LILRA3** is not on the GRCh38 primary assembly at all, so it is read twice:
as MAPQ-0 depth over the four alt contigs that carry it, and independently at the
deletion junction, where a LILRA3-bearing chromosome's reads soft-clip and a
deleted one's cross cleanly. The two share no failure mode — the depth assay dies
if the CRAM's reference lacks the alt contigs, and the junction assay is weakened
if the aligner *did* have them and moved the clipped reads there — so agreement
between them is worth more than either alone.

**The pair check.** LILRA6 copy number is also obtainable as the pooled
LILRA6+LILRB3 depth at MAPQ 0 minus LILRB3 from its unique window. It is a coarse
route — inverting the pooled ratio divides by the span weight and roughly doubles
the noise — so it is read for a systematic offset rather than per-sample
agreement. Disagreement means reads are moving between the paralogues.

Where this departs most from the predecessor is that copy number is **absolute
and per-sample**. `lilr-genotyper` fits the per-haploid depth unit across a whole
cohort, which makes ``cn_cohort`` a barrier in the DAG and means a cohort of
fewer than ~30 samples cannot be called at all. With WGS, λ₁ is measured from
each sample's own control loci, so one sample is enough. :func:`refine_cohort`
survives as an optional check on the scale, not as a prerequisite.
"""

from __future__ import annotations

import json
import math
import re
import statistics
from dataclasses import dataclass, field

from . import loci
from .coverage import CoverageModel
from .shell import require, run

# Copy-number ranges seen in pangenome data. Estimates are clamped to these, so
# an estimate far outside becomes a low-confidence call at the boundary rather
# than an impossible integer.
CN_RANGE = {
    "LILRA3": (0, 2),      # biallelic insertion/deletion
    "LILRA6": (0, 6),
    "LILRB3": (0, 4),
}

# An estimate this close to a half-integer boundary is a coin flip between the
# two neighbouring integers. Calls inside the band are still made — refusing
# would bias allele frequencies toward whichever class is easier to measure —
# but they are flagged, and the flag is what the manual override file is for.
AMBIGUOUS_BAND = 0.20

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


@dataclass
class CNCall:
    """One gene's copy number in one sample, and how much to believe it."""

    sample: str
    gene: str
    estimate: float | None = None        # continuous, before rounding
    copies: int | None = None
    confidence: float = 0.0              # 0 at a boundary, 1 at an integer
    method: str = ""
    status: str = "not_measured"         # measured | not_measured | failed
    support: dict = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)

    @property
    def ambiguous(self) -> bool:
        return self.status == "measured" and self.confidence < (1 - 2 * AMBIGUOUS_BAND)

    def as_row(self) -> dict:
        return {
            "sample": self.sample,
            "gene": self.gene,
            "copies": "" if self.copies is None else self.copies,
            "estimate": "" if self.estimate is None else round(self.estimate, 3),
            "confidence": round(self.confidence, 3),
            "method": self.method,
            "status": self.status,
            "ambiguous": self.ambiguous,
            "notes": ";".join(self.notes),
            # The raw counts the call was made from — depths, junction reads, the
            # pooled-pair figures. Carried because a copy number without its
            # evidence cannot be re-adjudicated later, and in this cluster the
            # calls that need re-adjudicating are exactly the plausible ones.
            "support": json.dumps(self.support, sort_keys=True),
        }


def integerise(estimate: float, gene: str) -> tuple[int, float]:
    """Round an estimate to a copy number, with a confidence.

    Confidence is the distance from the nearest half-integer boundary, scaled so
    that landing exactly on an integer gives 1.0 and landing on a boundary gives
    0.0. It measures only how cleanly this sample's estimate rounds — not whether
    λ₁ was right, whether the alignment was ALT-aware, or whether the gene was
    measurable at all. Those are separate fields precisely so that a confident
    number derived from a broken measurement cannot masquerade as a good call.
    """
    lo, hi = CN_RANGE.get(gene, (0, 6))
    nearest = round(estimate)
    copies = int(min(hi, max(lo, nearest)))
    confidence = 1.0 - 2.0 * abs(estimate - nearest)

    # The penalty applies only when clamping actually changed the answer. An
    # estimate of 2.047 at a gene capped at CN 2 rounds to 2 either way and is a
    # good call; treating "outside the range" as "outside by any amount" scored
    # four of five correct LILRA3 calls at zero confidence, which would flag them
    # as ambiguous and -- since the flag fires unevenly across copy-number
    # classes -- bias any allele frequency computed from the confident subset.
    if nearest != copies:
        confidence = 0.0
    return copies, max(0.0, confidence)


def _mean_depth(bam: str, intervals: list[tuple[str, int, int]], mapq: int,
                *, reference: str | None = None, samtools: str = "samtools",
                supplementary: bool = False,
                ) -> tuple[float, int] | None:
    """Mean per-base depth over a set of intervals, and the bases behind it.

    Args:
        supplementary: count supplementary alignments. False everywhere except
            LILRA3, where they are the entire signal — see
            :func:`_call_lilra3`.

    Returns None — never 0.0 — when the measurement could not be made. At
    LILRA3, zero is the *expected* answer for a deletion homozygote at ~24%
    allele frequency, so conflating "no reads" with "could not ask" would
    manufacture deletions out of failed queries.
    """
    require(samtools)
    bed = "".join(f"{c}\t{s}\t{e}\n" for c, s, e in intervals)
    # -G excludes a flag class. Built here rather than spliced in afterwards:
    # inserting it at the wrong index puts it before the `depth` subcommand,
    # which samtools rejects — and since a failed query returns None, the
    # symptom was every LILRA6 and LILRB3 call coming back "failed" rather than
    # anything pointing at a malformed command.
    exclude = [] if supplementary else ["-G", "0x800"]
    cmd = [samtools, "depth", "-a", "-Q", str(mapq), *exclude, "-b", "/dev/stdin"]
    if reference:
        cmd += ["--reference", reference]
    cmd.append(bam)
    try:
        out = run(cmd, text_input=bed).stdout
    except Exception:
        return None

    total = n = 0
    for line in out.splitlines():
        parts = line.split("\t")
        if len(parts) >= 3:
            total += int(parts[2])
            n += 1
    if n == 0:
        return None
    return total / n, n


# A clip has to be long enough to be sequence rather than a trimmed base or two,
# and the aligned part long enough to place it. Both at 10 bp, which also makes
# the two read classes comparable: a spanning read is required to hold 10 aligned
# bases either side of the breakpoint, a clipped one 10 aligned and 10 clipped.
MIN_CLIP = MIN_ANCHOR = 10
# The breakpoint's two ends are 5 bp apart, and where a read's clip is assigned
# within that depends on the aligner's tie-break, so allow a few bases either way.
JUNCTION_TOLERANCE = 5


def junction_counts(bam: str, *, reference: str | None = None,
                    samtools: str = "samtools", window: int = 300,
                    ) -> tuple[int, int] | None:
    """Reads clipped at the LILRA3 deletion junction, and reads spanning it.

    The primary assembly carries the deleted allele, so a chromosome that is also
    deleted matches the reference and its reads cross the junction cleanly, while
    one carrying LILRA3 diverges there and its reads soft-clip. The clipped
    fraction estimates half the copy number.

    The breakpoint has two ends, 5 bp apart across the microhomology: reads
    running into LILRA3 from the left flank end right-clipped at the first, reads
    coming back out of it start left-clipped at the second. Both are the bearing
    chromosome's signal, so both count as clipped.

    Returns (clipped, spanning), or None if the query failed.
    """
    require(samtools)
    chrom, left = loci.LILRA3_JUNCTION
    right = left + loci.LILRA3_JUNCTION_MICROHOMOLOGY
    region = f"{chrom}:{max(1, left - window)}-{right + window}"
    cmd = [samtools, "view", "-q", "1"]
    if reference:
        cmd += ["--reference", reference]
    cmd += [bam, region]
    try:
        out = run(cmd).stdout
    except Exception:
        return None

    clipped, spanning = tally_junction(out.splitlines(), left, right)
    if clipped + spanning == 0:
        return None
    return clipped, spanning


def tally_junction(sam_lines, left: int, right: int) -> tuple[int, int]:
    """Sort SAM records into breakpoint-clipped and breakpoint-spanning.

    Split out from :func:`junction_counts` so the classification can be tested
    on SAM text rather than needing a BAM and a samtools.
    """
    clipped = spanning = 0
    for line in sam_lines:
        fields = line.split("\t")
        if len(fields) < 6:
            continue
        start = int(fields[3]) - 1
        ops = CIGAR.findall(fields[5])
        if not ops:
            continue
        ref_len = sum(int(n) for n, op in ops if op in "MDN=X")
        end = start + ref_len
        if ref_len < MIN_ANCHOR:
            continue
        head, tail = ops[0], ops[-1]
        # Soft-clipping at the breakpoint is the signal; a clip elsewhere in the
        # read is ordinary adapter or quality trimming.
        if (tail[1] == "S" and int(tail[0]) >= MIN_CLIP
                and abs(end - left) <= JUNCTION_TOLERANCE):
            clipped += 1
        elif (head[1] == "S" and int(head[0]) >= MIN_CLIP
                and abs(start - right) <= JUNCTION_TOLERANCE):
            clipped += 1
        elif start < left - MIN_ANCHOR and end > right + MIN_ANCHOR:
            spanning += 1
    return clipped, spanning


def call_sample(sample: str, bam: str, model: CoverageModel, *,
                reference: str | None = None, samtools: str = "samtools",
                alt_depth_valid: bool = True,
                ) -> list[CNCall]:
    """Copy number for the three variable genes in one sample.

    Args:
        alt_depth_valid: whether this BAM contains every read that aligned to the
            LRC alt contigs. True for a slice taken straight out of a CRAM. False
            for a BAM built by realigning a regional extraction — see
            :func:`_call_lilra3` for what that costs and why the junction assay
            is used instead.
    """
    calls = [
        _call_unique_window(sample, "LILRA6", bam, model,
                            reference=reference, samtools=samtools),
        _call_unique_window(sample, "LILRB3", bam, model,
                            reference=reference, samtools=samtools),
        _call_lilra3(sample, bam, model, reference=reference, samtools=samtools,
                     alt_depth_valid=alt_depth_valid),
    ]
    _pair_check(calls, bam, model, reference=reference, samtools=samtools)
    return calls


def _call_unique_window(sample: str, gene: str, bam: str, model: CoverageModel,
                        *, reference: str | None, samtools: str) -> CNCall:
    call = CNCall(sample=sample, gene=gene, method="unique_window_q20")

    if model.lambda1 <= 0:
        call.status = "failed"
        call.notes.append("no coverage model")
        return call

    if not model.usable_mapq20:
        # The refusal that matters. A dead MAPQ-20 count over a live baseline is
        # a clean zero, and a cohort of those is a cohort of deletion
        # homozygotes reported with confidence.
        call.status = "not_measured"
        call.notes.append(
            f"MAPQ 20 is not usable in the LRC (q20_lrc={model.q20_lrc:.3f}, "
            f"verdict={model.alt_verdict}); reporting not measured rather than zero"
        )
        return call

    measured = _mean_depth(bam, loci.UNIQUE_WINDOWS[gene], loci.MAPQ_STRICT,
                           reference=reference, samtools=samtools)
    if measured is None:
        call.status = "failed"
        call.notes.append("depth query failed")
        return call

    depth, n_bases = measured
    call.estimate = depth / model.lambda1
    call.copies, call.confidence = integerise(call.estimate, gene)
    call.status = "measured"
    call.support = {"mean_depth": round(depth, 2), "n_bases": n_bases,
                    "lambda1": round(model.lambda1, 2)}
    return call


def _call_lilra3(sample: str, bam: str, model: CoverageModel, *,
                 reference: str | None, samtools: str,
                 alt_depth_valid: bool = True) -> CNCall:
    """LILRA3 by two independent routes, preferring depth and checking it.

    Depth is preferred because it is a direct measurement over 6.7 kb, where the
    junction assay rests on clipping behaviour at a breakpoint that sits inside
    an Alu. But the depth route needs the alt contigs in the CRAM's reference,
    and where they are absent the junction is all there is — so the fallback is
    real, not decorative.

    ``alt_depth_valid=False`` is the second reason that route can be unavailable,
    and it is not visible in the BAM the way a missing contig is. The MAPQ-0
    alt-contig depth counts supplementary records, and on HG00138's CRAM slice
    135 of the 964 reads with an alt-contig record have **no primary record in
    the LRC at all** — their primaries are scattered over chr2, chr3, chrX and
    the rest of the genome, repeat-derived reads with a supplementary hit on the
    LILRA3 contigs. A BAM built by extracting the LRC and realigning cannot
    contain them: the FASTQ step drops supplementary records, because emitting
    one would write a read twice, so a read whose *only* slice record is
    supplementary disappears.

    That is a 21% loss of alt-contig records on HG00138 (1,825 -> 1,448) and it
    is systematic: across 88 samples the realigned estimate is a median 0.404
    copies below the CRAM-as-is one and never above it, turning 19 true CN 2
    calls into CN 1. The depth route's calibration includes that repeat-derived
    component; a regional extraction excludes it; the same threshold cannot
    serve both. The junction assay is unaffected — it reads clipping at
    chr19:54,297,005, which is inside any LRC extraction — and reproduces
    itself across the two (HG00138: 47 clipped/0 spanning as-is, 45/0
    realigned, estimate 2.00 both ways). So it is used outright rather than as a
    cross-check.
    """
    call = CNCall(sample=sample, gene="LILRA3", method="alt_depth_q0")

    junction = junction_counts(bam, reference=reference, samtools=samtools)
    junction_estimate = None
    if junction is not None:
        clipped, spanning = junction
        # Over the 101-donor overlap this reads 0.00 at truth CN 0, a median 1.12
        # at CN 1 and exactly 2.00 at CN 2 — the 12% at CN 1 is the bearing
        # chromosome offering two breakpoints' worth of clipped reads against the
        # deleted one's single spanning window. Left uncorrected: it is well
        # inside the rounding band, and a fitted fudge factor on a cross-check
        # would couple it to the route it exists to be independent of.
        junction_estimate = 2.0 * clipped / (clipped + spanning)
        call.support["junction_clipped"] = clipped
        call.support["junction_spanning"] = spanning
        call.support["junction_estimate"] = round(junction_estimate, 3)

    # MAPQ 0, and normalised on a MAPQ-0 baseline: a ratio has to divide like by
    # like, and dividing a MAPQ-0 count by a MAPQ-20 baseline would inflate it by
    # whatever fraction of the baseline's reads are repeat-derived.
    lrc_q0 = [c.mean_q0 for c in model.controls if c.inside_placement and c.mean_q0 > 0]
    # supplementary=True: `bwa mem -Y` emits the ALT-contig hit of an ALT-aware
    # alignment as a supplementary record, and since LILRA3 is absent from the
    # primary assembly that is where essentially all of its evidence sits. On
    # HG00099, 481 of 547 records over one alt interval are supplementary;
    # excluding them made a two-copy donor read as a deletion homozygote.
    measured = None
    if alt_depth_valid:
        measured = _mean_depth(bam, loci.LILRA3_ALT, loci.MAPQ_ANY,
                               reference=reference, samtools=samtools,
                               supplementary=True)

    if measured is not None and lrc_q0:
        depth, n_bases = measured
        # Sum across the four intervals, divide by ONE interval's length. A read
        # from a LILRA3-bearing chromosome lands on exactly one of the four --
        # measured on HG00099, two of the contigs share zero read names out of
        # 518 and 464 -- so the four counts partition the evidence rather than
        # replicating it. `_mean_depth` averages over all four intervals, so the
        # per-interval mean has to be multiplied back up by their number.
        depth = depth * len(loci.LILRA3_ALT)
        baseline_haploid = statistics.median(lrc_q0) / 2.0
        call.estimate = depth / baseline_haploid if baseline_haploid > 0 else None
        call.support["summed_depth"] = round(depth, 2)
        call.support["n_bases"] = n_bases
        call.support["n_alt_intervals"] = len(loci.LILRA3_ALT)

    if call.estimate is None and junction_estimate is not None:
        call.estimate = junction_estimate
        call.method = "junction"
        call.notes.append(
            "called from the deletion junction alone; the alt-contig depth route "
            + ("is not valid on a BAM built from a regional extraction, which "
               "cannot hold the genome-wide reads whose supplementary alignments "
               "that route counts"
               if not alt_depth_valid else
               "is unavailable (the CRAM's reference has no LRC alt contigs)")
        )

    if call.estimate is None:
        call.status = "failed"
        call.notes.append("neither the alt-contig depth nor the junction could be read")
        return call

    call.copies, call.confidence = integerise(call.estimate, "LILRA3")
    call.status = "measured"

    if junction_estimate is not None and call.method != "junction":
        disagreement = abs(junction_estimate - call.estimate)
        call.support["junction_disagreement"] = round(disagreement, 3)
        if disagreement > 0.75:
            call.notes.append(
                f"depth says {call.estimate:.2f} copies and the junction says "
                f"{junction_estimate:.2f}; these assays share no failure mode, so "
                "treat this call as unresolved"
            )
            call.confidence = min(call.confidence, 0.3)
    return call


def _pair_check(calls: list[CNCall], bam: str, model: CoverageModel, *,
                reference: str | None, samtools: str) -> None:
    """Cross-check LILRA6 against the pooled LILRA6+LILRB3 depth.

    Coarse by construction: inverting the pooled ratio divides by the span weight
    and roughly doubles the noise. Read it for a systematic offset across a
    cohort, not for per-sample agreement — which is why disagreement adds a note
    and does not change the call.
    """
    by_gene = {c.gene: c for c in calls}
    a6, b3 = by_gene.get("LILRA6"), by_gene.get("LILRB3")
    if not a6 or not b3 or a6.status != "measured" or b3.status != "measured":
        return
    if model.lambda1 <= 0:
        return

    gene_bodies = []
    spans = []
    for gene in ("LILRA6", "LILRB3"):
        chrom, start, end = loci.gene_span(gene)
        gene_bodies.append((chrom, start, end))
        spans.append(end - start)

    measured = _mean_depth(bam, gene_bodies, loci.MAPQ_ANY,
                           reference=reference, samtools=samtools)
    if measured is None:
        return
    depth, _ = measured

    # The span weight. Reads from both genes multi-map freely across both gene
    # bodies, so the pair's whole output spreads over the *sum* of the two spans
    # while one copy of one gene contributes coverage over one span. Mean depth
    # over the union is therefore the pair's copy number divided by two, not the
    # pair's copy number — and dropping the weight makes a 2+2 sample read as
    # "pooled = 2", from which the implied LILRA6 copy number comes out at zero
    # for everybody.
    total_span = sum(spans)
    mean_span = total_span / len(spans)
    pooled = depth * total_span / (model.lambda1 * mean_span)
    implied_a6 = pooled - (b3.copies or 0)
    a6.support["pooled_pair_estimate"] = round(pooled, 2)
    a6.support["pair_implied_lilra6"] = round(implied_a6, 2)
    if a6.copies is not None and abs(implied_a6 - a6.copies) > 1.5:
        a6.notes.append(
            f"the pooled pair implies ~{implied_a6:.1f} LILRA6 copies against "
            f"{a6.copies} from the unique window; reads may be moving between "
            "the paralogues"
        )


def refine_cohort(calls: list[CNCall]) -> dict:
    """Check the cohort's scale, and report rather than silently rescale.

    With an absolute baseline the per-copy unit should already be 1.0. If the
    cohort's estimates cluster around integers offset by a constant factor, that
    is a systematic error in λ₁ — a GC correction that is off, or control loci
    that are not at copy number 2 in this population — and it should be
    investigated, not divided out.

    This is deliberately not the predecessor's ``cn_cohort``. That fits the unit
    across the cohort and *uses* it, which makes the fit load-bearing and the
    cohort a barrier. Here the fit is a diagnostic and nothing downstream depends
    on it.
    """
    out: dict[str, dict] = {}
    by_gene: dict[str, list[float]] = {}
    for c in calls:
        if c.status == "measured" and c.estimate is not None:
            by_gene.setdefault(c.gene, []).append(c.estimate)

    for gene, estimates in by_gene.items():
        nonzero = [e for e in estimates if e > 0.25]
        if len(nonzero) < 20:
            out[gene] = {"n": len(estimates), "unit": None,
                         "note": "too few samples to check the scale"}
            continue
        if not _spans_enough_classes(nonzero):
            out[gene] = {
                "n": len(estimates), "unit": None,
                "median_estimate": round(statistics.median(nonzero), 3),
                "note": "copy number is too uniform in this cohort to fit a "
                        "unit; the spacing between classes is what identifies it",
            }
            continue
        unit = _fit_unit(nonzero)
        out[gene] = {
            "n": len(estimates),
            "unit": round(unit, 4),
            "median_estimate": round(statistics.median(nonzero), 3),
            "ambiguous_fraction": round(
                sum(1 for c in calls
                    if c.gene == gene and c.ambiguous) / max(len(estimates), 1), 4),
            "note": ("scale looks right" if abs(unit - 1.0) < 0.08 else
                     f"estimates cluster on a unit of {unit:.3f}, not 1.0 — "
                     "lambda1 is systematically off by that factor"),
        }
    return out


# A copy-number class has to carry this many samples before it counts as one of
# the two the unit fit is measured between.
MIN_CLASS_SUPPORT = 5


def _spans_enough_classes(estimates: list[float]) -> bool:
    """Whether the cohort has two populated copy-number classes to fit between.

    The unit is the *spacing* between classes, so a cohort sitting on one class
    does not constrain it: any unit that divides that one value near-integrally
    fits about as well as any other. LILRB3 in the 101-donor overlap is 99 donors
    at CN 2 and 2 at CN 1, and the grid duly reported a unit of 0.700 — cost
    9.817 against 9.973 for the correct 1.03, a 1.6% margin — with the note
    "lambda1 is systematically off by that factor", on the one gene that scored
    100% against truth. A diagnostic that cries wolf on a perfect result is worse
    than no diagnostic, so declining to fit is the honest answer here.
    """
    counts: dict[int, int] = {}
    for e in estimates:
        counts[round(e)] = counts.get(round(e), 0) + 1
    return sum(1 for n in counts.values() if n >= MIN_CLASS_SUPPORT) >= 2


def _fit_unit(estimates: list[float]) -> float:
    """The spacing that the cohort's estimates actually cluster on.

    Minimises total rounding error over a grid. A grid rather than an optimiser
    keeps it dependency-free and the objective is not smooth enough for the
    derivative-based alternative to be trustworthy anyway.
    """
    best_unit, best_cost = 1.0, math.inf
    for i in range(700, 1401):
        unit = i / 1000.0
        cost = sum(abs(e - round(e / unit) * unit) for e in estimates)
        if cost < best_cost:
            best_unit, best_cost = unit, cost
    return best_unit
