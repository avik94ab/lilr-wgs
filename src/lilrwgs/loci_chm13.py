"""T2T-CHM13v2.0 coordinates for the leukocyte receptor complex.

The same intervals as :mod:`lilrwgs.loci`, in the assembly that actually contains
LILRA3. Every value here was derived by ``scripts/derive_loci_chm13.py`` from the
CHM13 sequence and annotation, never lifted over: liftover is least trustworthy
in exactly this kind of region, and a coordinate that lands one paralogue to the
left produces copy numbers rather than errors.

**LILRA3 is on the primary assembly here, and that is the whole point.** GRCh38's
chr19 carries the common ~6.7 kb deletion, so LILRA3 lives only on alt contigs
and has to be read as MAPQ-0 depth across four near-identical haplotypes while
counting supplementary records — a route that cannot survive a regional
extraction, because the reads it counts have primaries scattered genome-wide
(PLAN.md §12). CHM13 carries the insertion allele, so LILRA3 is ordinary
single-copy sequence measurable at MAPQ 20 like any other gene.

That was measured, not assumed. The 7,126 bp calling reference
``resources/bundle/references/LILRA3_named.fa`` aligns here at 7,126/7,126
identity, reverse-complemented, MAPQ 60, and every other hit in the region is a
partial paralogue at 59–87%. Note the CHM13 annotation does **not** list
LILRA3: ``chm13v2.0_RefSeq_Liftoff_v5.1.gff3`` is lifted from the GRCh38 primary
assembly, which has no LILRA3 to lift. Its absence there is a fact about the
annotation's provenance and says nothing about the sequence.

**There are no ALT contigs.** CHM13 is a single complete haplotype, so the whole
apparatus that exists to survive GRCh38's nine LRC alt haplotypes placed at one
primary interval — the MAPQ collapse, :data:`lilrwgs.loci.ALT_PLACEMENT`, the
``.alt`` file that ``bwa`` reads without ever mentioning — has no counterpart and
needs none. :data:`OUTSIDE_LOCI` survives, but as an ordinary sanity check on
whether the region is being sampled normally rather than as an ALT-awareness
test.

**The cost of a single haplotype.** GRCh38-plus-alts represents several LRC
haplotypes; CHM13 represents one. In a region this polymorphic that may recruit
divergent haplotypes less well, and it is the reason this module is offered
beside :mod:`lilrwgs.loci` rather than replacing it. Compare the two on the same
samples before believing either.
"""

from __future__ import annotations

from .loci import (  # the assembly-independent parts
    MAPQ_ANY,
    MAPQ_STRICT,
    Control,
    Gene,
    as_region,
    merge_intervals,
)

__all__ = [
    "ASSEMBLY", "CHROM", "CHR19_LENGTH", "GENES", "GENE_BODIES", "LRC_SLICE",
    "LILRA3_SPAN", "UNIQUE_WINDOWS", "CONTROL_LOCI", "OUTSIDE_LOCI",
    "ALL_CONTROLS", "MAPQ_STRICT", "MAPQ_ANY", "FLANK", "Control", "Gene",
    "as_region", "merge_intervals", "gene_span", "slice_intervals",
    "slice_regions", "extraction_regions",
]

ASSEMBLY = "CHM13v2.0"
CHROM = "chr19"
# Used to tell this assembly from GRCh38, whose chr19 is 58,617,616.
CHR19_LENGTH = 61_707_364

# 5' to 3' along the cluster, the same order as GRCh38 — including LILRA3, which
# here is a gene like any other rather than a special case.
GENES = [
    "LILRB3", "LILRA6", "LILRB5", "LILRB2", "LILRA3", "LILRA5",
    "LILRA4", "LILRA2", "LILRA1", "LILRB1", "LILRB4",
]

# RefSeq Liftoff v5.1 gene spans, except LILRA3, which that file does not carry
# and which was located by aligning the calling reference (identity 1.0000).
GENE_BODIES: dict[str, Gene] = {
    g.name: g
    for g in [
        Gene("LILRB3", CHROM, 57_296_665, 57_303_395, "-"),
        Gene("LILRA6", CHROM, 57_316_984, 57_323_158, "-"),
        Gene("LILRB5", CHROM, 57_331_817, 57_340_006, "-"),
        Gene("LILRB2", CHROM, 57_356_203, 57_363_500, "-"),
        Gene("LILRA3", CHROM, 57_377_083, 57_384_209, "-",
             "present on the primary assembly, unlike GRCh38; located by "
             "aligning LILRA3_named.fa at identity 1.0000"),
        Gene("LILRA5", CHROM, 57_396_883, 57_402_937, "-"),
        Gene("LILRA4", CHROM, 57_422_955, 57_428_933, "-"),
        Gene("LILRA2", CHROM, 57_666_504, 57_683_802, "+"),
        Gene("LILRA1", CHROM, 57_687_146, 57_695_896, "+"),
        Gene("LILRB1", CHROM, 57_709_826, 57_731_530, "+"),
        Gene("LILRB4", CHROM, 57_756_497, 57_762_244, "+"),
    ]
}

# LILRA3's span, for symmetry with lilrwgs.loci.LILRA3_SPAN. There it is one
# interval's length shared across four alt contigs; here it is simply the gene.
LILRA3_SPAN = 7_126

# ---------------------------------------------------------------------------
# Extraction
# ---------------------------------------------------------------------------

# Bounded below by LILRB3 with ~67 kb of flank and above by LILRB4 with ~31 kb,
# the same geometry as the GRCh38 slice. The upper bound stops short of KIR3DL3
# at chr19:57,817,388 — measured in this assembly, not extrapolated — leaving the
# same ~24 kb margin GRCh38 leaves. The KIR cluster is the most copy-number
# variable region in the genome and pulling it in would multiply the extracted
# readset for nothing.
LRC_SLICE = (CHROM, 57_230_000, 57_793_000)

# ---------------------------------------------------------------------------
# Paralogue-unique windows
# ---------------------------------------------------------------------------

# Transferred from lilrwgs.loci.UNIQUE_WINDOWS by aligning their sequence, not
# re-derived. Re-running the documented criterion — a 100-mer is ambiguous if
# another within 3 mismatches exists anywhere in the LRC on either strand —
# against GRCh38 reproduces these windows' *starts* exactly but runs longer and
# finds extra windows, 4,003 bp against 2,900 at LILRA6. So the published
# criterion is looser than whatever produced the windows LILRA6's accuracy was
# validated on, and transferring the validated sequence keeps the definition that
# earned that accuracy while changing only the assembly it is expressed in.
#
# All four transferred at >=0.9975 identity with their lengths preserved exactly,
# and the independent 100-mer audit finds them still unique here: 0.0%, 0.0%,
# 0.4% and 0.0% ambiguous. Both numbers are in resources/chm13_loci.tsv.
UNIQUE_WINDOWS = {
    "LILRA6": [(CHROM, 57_316_984, 57_319_070),    # 2,086 bp, identity 1.0000
               (CHROM, 57_319_444, 57_320_258)],   # 814 bp,   identity 0.9975
    "LILRB3": [(CHROM, 57_296_847, 57_297_439),    # 592 bp,   identity 1.0000
               (CHROM, 57_298_639, 57_299_928)],   # 1,289 bp, identity 0.9984
    # LILRA3 is the exception, and it has to be derived rather than transferred:
    # there is no GRCh38 window to move, because GRCh38 has no LILRA3.
    #
    # It needs no window in the usual sense. LILRA6 and LILRB3 are 47% and 60%
    # ambiguous by the 100-mer criterion, which is why only their 3' ends are
    # measurable; LILRA3 is **9.6%**, and the ambiguity is scattered — 677
    # ambiguous 100-mer starts in 34 runs, the longest spanning ~182 bp, with no
    # contiguous dead zone. Every base of the gene is covered by some
    # unambiguous 100-mer, and since a read is 150 bp and the criterion's window
    # is 100, a read over any base can carry unique sequence.
    #
    # That is biology rather than luck: LILRA3 is the soluble family member,
    # lacking the transmembrane and cytoplasmic domains, and is not a recent
    # duplicate of a neighbour the way LILRA6 and LILRB3 are of each other.
    #
    # So the window is the gene. This is the measurement GRCh38 cannot make at
    # all — there, LILRA3 is MAPQ-0 depth over four alt haplotypes.
    "LILRA3": [(CHROM, 57_377_083, 57_384_209)],   # 7,126 bp, 9.6% ambiguous
}

# ---------------------------------------------------------------------------
# Controls
# ---------------------------------------------------------------------------

# The same housekeeping loci as GRCh38, at their CHM13 spans. `ambiguous` is
# carried over from the GRCh38 measurement and marked as such: it is documentary
# — nothing in this package branches on it — and re-measuring it here would say
# something about CHM13's duplication structure rather than about the pipeline.
CONTROL_LOCI = [
    Control("PRPF31", CHROM, 57_194_450, 57_210_410, 0.017, True),
    Control("TMC4",   CHROM, 57_240_373, 57_253_360, 0.000, True),
    Control("MBOAT7", CHROM, 57_253_614, 57_269_968, 0.004, True),
    Control("TTYH1",  CHROM, 57_508_740, 57_530_183, 0.001, True),
]

# Controls further from the cluster. On GRCh38 these answer "was the alignment
# ALT-aware?", a question CHM13 cannot pose because it has no ALT contigs. They
# are kept because the other half of their job still applies: if MAPQ-20
# retention or depth differs between the LRC's neighbourhood and sequence a
# megabase away, something about this sample's sampling is wrong.
OUTSIDE_LOCI = [
    Control("PRKCG",  CHROM, 56_960_798, 56_987_358, 0.004, False),
    Control("CACNG6", CHROM, 57_069_638, 57_091_148, 0.000, False),
    Control("PPP6R1", CHROM, 58_322_177, 58_352_962, 0.002, False),
]

ALL_CONTROLS = CONTROL_LOCI + OUTSIDE_LOCI

# Flank added around every fetched interval; see lilrwgs.loci.FLANK.
FLANK = 1_000


def gene_span(name: str) -> tuple[str, int, int]:
    """Primary-assembly interval of a gene body.

    Unlike :func:`lilrwgs.loci.gene_span`, this never raises for LILRA3 — the
    entire reason this module exists is that CHM13 carries it.
    """
    g = GENE_BODIES[name]
    return (g.chrom, g.start, g.end)


def slice_intervals(include_alts: bool = True) -> list[tuple[str, int, int]]:
    """One pass per sample: everything any stage will need, merged.

    ``include_alts`` is accepted and ignored, so this module is a drop-in for
    :mod:`lilrwgs.loci`. CHM13 has no alt contigs and LILRA3 is inside
    :data:`LRC_SLICE`, so there is nothing for the flag to include.
    """
    intervals = [LRC_SLICE]
    intervals.extend((c.chrom, c.start - FLANK, c.end + FLANK)
                     for c in ALL_CONTROLS)
    return merge_intervals([(c, max(0, s), e) for c, s, e in intervals])


def slice_regions(include_alts: bool = True) -> list[str]:
    """:func:`slice_intervals` as samtools region strings."""
    return [as_region(*iv) for iv in slice_intervals(include_alts)]


def extraction_regions(include_alts: bool = True) -> list[str]:
    """Everywhere a LILR read could have been aligned in this assembly.

    Just the cluster: there are no alt contigs to add, and LILRA3 is in it.
    """
    return [as_region(*LRC_SLICE)]
