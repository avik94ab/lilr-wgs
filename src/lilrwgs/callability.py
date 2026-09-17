"""Per-position evidence, and the track that records it.

This is where :mod:`lilrwgs.depth_model` meets an actual BAM. For one
(sample, gene) it walks the alignment, measures three things at every position —
depth, the fraction of reads clearing the MAPQ floor, and the fraction that were
kept in more than one gene by read assignment — and asks the model whether that
position can be called.

The output is a *track*, not a mask. The predecessor wrote positions below its
threshold as ``N`` in the consensus and kept no record, so an N in its output
could equally mean no reads, reads that were ambiguous, or a real deletion.
Those are different findings and distinguishing them is most of the diagnostic
value: a gene that is 40% N because of low depth needs more coverage, and one
that is 40% N because of paralogue ambiguity will never be fixed by more reads.
"""

from __future__ import annotations

import csv
import os
from dataclasses import dataclass
from pathlib import Path

from .assign import SHARED_TAG
from .depth_model import Callability, PositionCall, call_position, summarise
from .loci import MAPQ_STRICT


@dataclass
class PositionEvidence:
    """What the reads say at one reference position."""

    pos: int                 # 0-based
    depth: int = 0
    n_pass_mapq: int = 0
    n_shared: int = 0

    @property
    def mapq_fraction(self) -> float:
        return self.n_pass_mapq / self.depth if self.depth else 0.0

    @property
    def shared_fraction(self) -> float:
        return self.n_shared / self.depth if self.depth else 0.0


def gather_evidence(bam: str | os.PathLike, ref_length: int, *,
                    shared_names: dict[str, int] | None = None,
                    mapq_floor: int = MAPQ_STRICT,
                    min_baseq: int = 13) -> list[PositionEvidence]:
    """Walk a per-locus BAM and count, per position, what supports it.

    Args:
        shared_names: read names kept across a shared group, from
            :func:`lilrwgs.assign.read_shared_names`. Reads are also checked for
            the ``ZS`` tag, so this works whether the BAM came straight from
            assignment or was realigned from FASTQ — the tag does not survive the
            FASTQ hop and the table is how the fact crosses it.
    """
    import pysam

    shared_names = shared_names or {}
    evidence = [PositionEvidence(pos=i) for i in range(ref_length)]

    with pysam.AlignmentFile(str(bam), "rb") as handle:
        for column in handle.pileup(stepper="nofilter", truncate=False,
                                    min_base_quality=min_baseq,
                                    max_depth=100_000):
            pos = column.reference_pos
            if not 0 <= pos < ref_length:
                continue
            slot = evidence[pos]
            for read in column.pileups:
                if read.is_del or read.is_refskip:
                    continue
                aln = read.alignment
                if aln.is_secondary or aln.is_supplementary:
                    continue
                slot.depth += 1
                if aln.mapping_quality >= mapq_floor:
                    slot.n_pass_mapq += 1
                if aln.has_tag(SHARED_TAG) or aln.query_name in shared_names:
                    slot.n_shared += 1
    return evidence


def classify(evidence: list[PositionEvidence], *, copies: int, lambda1: float,
             dispersion: float, paralog_copies: int = 0,
             alpha: float = 0.005, min_mapq_fraction: float = 0.5,
             gc_by_pos: list[float] | None = None,
             gc_lookup=None) -> list[PositionCall]:
    """Run the depth model over a gene's positions.

    Args:
        paralog_copies: copy number of the gene this one shares a block with, or
            0 if it is separable. Positions whose reads are mostly shared
            legitimately carry both genes' coverage, and without this the
            two-sided gate rejects them as over-covered — which would discard
            most of the LILRA6 coding sequence in every sample.
        gc_lookup: optional ``gc -> lambda`` callable, normally
            :meth:`lilrwgs.coverage.CoverageModel.lambda_at`. Applied per
            position so a GC-extreme exon is judged against what its own GC
            predicts rather than against the gene's average.
    """
    calls: list[PositionCall] = []
    for i, e in enumerate(evidence):
        lam = lambda1
        if gc_lookup is not None and gc_by_pos is not None and i < len(gc_by_pos):
            lam = gc_lookup(gc_by_pos[i])
        calls.append(call_position(
            e.depth, copies=copies, lambda1=lam, dispersion=dispersion,
            alpha=alpha, shared_fraction=e.shared_fraction,
            paralog_copies=paralog_copies, mapq_fraction=e.mapq_fraction,
            min_mapq_fraction=min_mapq_fraction,
        ))
    return calls


def write_track(path: str | os.PathLike, contig: str,
                evidence: list[PositionEvidence],
                calls: list[PositionCall]) -> dict:
    """Write the per-position track and return its summary.

    One row per position, with the numbers that produced the verdict, so a
    disputed call can be argued with rather than merely disbelieved.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["contig", "pos", "depth", "mapq_fraction", "shared_fraction",
                    "effective_copies", "expected", "floor", "ceiling", "status"])
        for e, c in zip(evidence, calls):
            w.writerow([contig, e.pos + 1, e.depth,
                        round(e.mapq_fraction, 3), round(e.shared_fraction, 3),
                        round(c.effective_copies, 2), round(c.expected, 1),
                        c.floor, c.ceiling, c.status.value])
    return summarise(calls)


def callable_bed(path: str | os.PathLike, contig: str,
                 calls: list[PositionCall]) -> int:
    """Write callable positions as a BED of merged runs.

    This is what constrains variant calling. Where the predecessor passed
    ``bcftools view -i 'FMT/DP>=20'`` — one constant for a whole gene — the
    equivalent here is a region file whose bounds vary per position with copy
    number, local GC and whether the block is shared. Feeding it to
    ``bcftools view -T`` applies the whole model without any of it having to be
    expressed as a filter expression.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    runs: list[tuple[int, int]] = []
    start = None
    for i, c in enumerate(calls):
        if c.status is Callability.OK:
            if start is None:
                start = i
        elif start is not None:
            runs.append((start, i))
            start = None
    if start is not None:
        runs.append((start, len(calls)))

    with path.open("w") as fh:
        for s, e in runs:
            fh.write(f"{contig}\t{s}\t{e}\n")
    return sum(e - s for s, e in runs)


def mask_sequence(seq: str, calls: list[PositionCall]) -> str:
    """Mask uncallable positions to N.

    Deliberately paired with :func:`write_track`: the N tells a consumer the base
    is unknown, and the track tells them why. Emitting one without the other is
    what makes an N-heavy consensus uninterpretable.
    """
    if not seq:
        return seq
    out = list(seq)
    for i, c in enumerate(calls):
        if i >= len(out):
            break
        if c.status is not Callability.OK:
            out[i] = "N"
    return "".join(out)


def gc_by_position(ref_seq: str, window: int = 100) -> list[float]:
    """Local GC fraction at every position of a reference.

    Windowed at the read length, which is the scale at which GC affects whether
    a fragment was amplified and sequenced at all.
    """
    n = len(ref_seq)
    seq = ref_seq.upper()
    half = window // 2
    out: list[float] = []
    for i in range(n):
        lo, hi = max(0, i - half), min(n, i + half)
        chunk = seq[lo:hi]
        acgt = sum(chunk.count(b) for b in "ACGT")
        out.append((chunk.count("G") + chunk.count("C")) / acgt if acgt else 0.0)
    return out
