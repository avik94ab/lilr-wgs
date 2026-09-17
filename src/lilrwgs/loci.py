"""GRCh38 coordinates for the leukocyte receptor complex.

Every interval here is half-open and 0-based unless the name says otherwise, and
every one is traceable to a source rather than to memory:

* Gene bodies come from UCSC ``ncbiRefSeqCurated`` on hg38, taking the union of
  transcript spans per gene symbol (``scripts/derive_loci.py`` regenerates this
  file and diffs it against what is here).
* The paralogue-unique windows, the LILRA3 alt-contig intervals, the deletion
  junction and the control loci come from ``lilrCN_aou`` , which derived them by
  asking, for every 100-mer of each gene, whether a twin within 3 mismatches
  exists anywhere in chr19:53.9-55.1 Mb on either strand.

Two facts about this region drive most of the design and are worth stating where
the coordinates live:

**LILRA3 is not on the primary assembly.** The GRCh38 reference chromosome
carries the common ~6.7 kb deletion, so LILRA3 is annotated only on LRC alt
contigs. Depth over LILRA3 is depth over those contigs, and because they are
near-identical to each other those reads multi-map and exist only at MAPQ 0.
That is a deliberate exception to the MAPQ floor used everywhere else.

**Nine alt contigs are all placed at one primary interval** that contains every
window here, targets and controls alike. Given the ``.alt`` file, bwa computes
MAPQ over primary hits alone and a read in an LRC window keeps it; without it,
the same read has nine more equally good placements and comes back MAPQ 0. So
every MAPQ-20 number in this package is conditional on the alignment having been
ALT-aware, which is why :data:`OUTSIDE_LOCI` exists and why nothing trusts it
without measuring it.
"""

from __future__ import annotations

from dataclasses import dataclass, field

ASSEMBLY = "GRCh38"
CHROM = "chr19"

# The 11 target genes, 5' to 3' along the cluster as it is laid out on chr19.
GENES = [
    "LILRB3", "LILRA6", "LILRB5", "LILRB2", "LILRA3", "LILRA5",
    "LILRA4", "LILRA2", "LILRA1", "LILRB1", "LILRB4",
]


@dataclass(frozen=True)
class Gene:
    """A LILR gene body on the primary assembly.

    ``start``/``end`` are None for LILRA3, which the primary assembly does not
    carry; :data:`LILRA3_ALT` holds where it actually lives.
    """

    name: str
    chrom: str
    start: int | None
    end: int | None
    strand: str
    note: str = ""

    @property
    def on_primary(self) -> bool:
        return self.start is not None


# UCSC ncbiRefSeqCurated, hg38, union of curated transcript spans per symbol.
GENE_BODIES: dict[str, Gene] = {
    g.name: g
    for g in [
        Gene("LILRB3", CHROM, 54_216_277, 54_223_007, "-"),
        Gene("LILRA6", CHROM, 54_236_589, 54_242_790, "-"),
        Gene("LILRB5", CHROM, 54_249_420, 54_257_273, "-"),
        Gene("LILRB2", CHROM, 54_273_811, 54_281_110, "-"),
        Gene("LILRA3", CHROM, None, None, "-",
             "deleted from the primary assembly; see LILRA3_ALT"),
        Gene("LILRA5", CHROM, 54_307_069, 54_313_166, "-"),
        Gene("LILRA4", CHROM, 54_333_184, 54_339_162, "-"),
        Gene("LILRA2", CHROM, 54_572_987, 54_590_287, "+"),
        Gene("LILRA1", CHROM, 54_593_631, 54_602_381, "+"),
        Gene("LILRB1", CHROM, 54_616_321, 54_638_022, "+"),
        Gene("LILRB4", CHROM, 54_662_984, 54_668_718, "+"),
    ]
}

# ---------------------------------------------------------------------------
# Extraction
# ---------------------------------------------------------------------------

# The slice taken out of each CRAM. Bounded below by LILRB3 and above by LILRB4
# with ~65 kb of flank either side, which is far more than a 450 bp insert needs
# to keep pairs intact at the edges. The upper bound deliberately stops short of
# KIR3DL3 at chr19:54,724,441: the KIR cluster is the most copy-number-variable
# region in the genome and pulling it in would multiply the extracted readset for
# nothing.
LRC_SLICE = (CHROM, 54_150_000, 54_700_000)

# The nine LRC alt haplotypes, all placed by UCSC at one primary interval. Reads
# from a LILRA3-bearing chromosome have their only primary alignment on one of
# these, so extraction that skipped them would silently drop LILRA3 entirely.
ALT_CONTIGS = [
    "chr19_GL949746v1_alt", "chr19_GL949747v2_alt", "chr19_GL949748v2_alt",
    "chr19_GL949749v2_alt", "chr19_GL949750v2_alt", "chr19_GL949751v2_alt",
    "chr19_GL949752v1_alt", "chr19_GL949753v2_alt", "chr19_KI270938v1_alt",
]
ALT_PLACEMENT = (CHROM, 54_025_633, 55_084_318)

# ---------------------------------------------------------------------------
# LILRA3: the deletion
# ---------------------------------------------------------------------------

# The four alt contigs that carry LILRA3, and the interval on each. Anchoring
# each to primary chr19 on unique 41-mers puts the breakpoint at the same base on
# all four, at 0.98 identity across the junction neighbourhood against 0.63 for
# either flank alone.
LILRA3_ALT = [
    ("chr19_GL949746v1_alt", 271_714, 278_472),
    ("chr19_GL949747v2_alt", 271_595, 278_354),
    ("chr19_GL949753v2_alt", 271_970, 278_728),
    ("chr19_KI270938v1_alt", 271_943, 278_701),
]
# One copy, not the sum of four: a read from a LILRA3-bearing chromosome has
# exactly one primary alignment spread across four near-identical intervals, so
# the counts are summed and divided by a single interval's length.
LILRA3_SPAN = 6_758

# Primary is the deleted allele, so a deleted chromosome's reads cross this base
# cleanly and a LILRA3-bearing one's soft-clip there. clipped/(clipped+spanning)
# estimates half the copy number — an assay with no failure mode in common with
# the depth route.
LILRA3_JUNCTION = (CHROM, 54_296_977)

# ~980 bp of LILRA3's 3' end survives the deletion and sits on primary just left
# of the junction, at two copies in everyone. A QC probe, not a target: if it
# does not read flat near 1.0, the region is not being sampled as assumed.
LILRA3_RETAINED = [(CHROM, 54_295_997, 54_296_977)]

# ---------------------------------------------------------------------------
# Paralogue-unique windows
# ---------------------------------------------------------------------------

# LILRA6 and LILRB3 are ~97% identical across the 5' exons encoding the
# extracellular Ig domains. The unique sequence is all at the 3' end, encoding
# the transmembrane and cytoplasmic regions, where an activating receptor and an
# inhibitory one genuinely differ. A whole-gene window at MAPQ 20 would measure
# mostly nothing; these are where a MAPQ floor means something.
UNIQUE_WINDOWS = {
    "LILRA6": [(CHROM, 54_236_589, 54_238_675), (CHROM, 54_239_049, 54_239_863)],
    "LILRB3": [(CHROM, 54_216_459, 54_217_051), (CHROM, 54_218_251, 54_219_540)],
}

# ---------------------------------------------------------------------------
# Controls
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Control:
    name: str
    chrom: str
    start: int
    end: int
    ambiguous: float       # fraction of 100-mers with a twin, primary only
    inside_placement: bool

    @property
    def span(self) -> int:
        return self.end - self.start


# Housekeeping loci within ~200 kb of the targets, none a known CNV, all measured
# clean by the same 100-mer criterion. LAIR1 sits closer and was the obvious
# fifth, but it is 28.4% ambiguous against LAIR2 and is deliberately not used.
CONTROL_LOCI = [
    Control("PRPF31", CHROM, 54_115_410, 54_131_719, 0.017, True),
    Control("TMC4",   CHROM, 54_160_095, 54_173_250, 0.000, True),
    Control("MBOAT7", CHROM, 54_173_412, 54_189_882, 0.004, True),
    Control("TTYH1",  CHROM, 54_415_219, 54_436_904, 0.001, True),
]

# Controls of the same kind placed *outside* the alt placement, so no alt contig
# carries a second copy of them however bwa was run. They answer a question the
# LRC controls cannot: was the alignment ALT-aware? With the alt contigs in the
# search space, 100% of LILRA6 and LILRB3 100-mers have a twin and the four LRC
# controls go from 0.0-1.7% ambiguous to 96-99%. These do not move.
OUTSIDE_LOCI = [
    Control("PRKCG",  CHROM, 53_882_196, 53_907_652, 0.004, False),
    Control("CACNG6", CHROM, 53_991_148, 54_012_666, 0.000, False),
    Control("PPP6R1", CHROM, 55_229_778, 55_259_017, 0.002, False),
]

ALL_CONTROLS = CONTROL_LOCI + OUTSIDE_LOCI

# MAPQ floors. 20 everywhere a floor can work; 0 where the sequence is
# multi-mapping by construction and a floor would measure nothing.
MAPQ_STRICT = 20
MAPQ_ANY = 0


def inside_alt_placement(chrom: str, start: int, end: int) -> bool:
    """Whether an interval overlaps the primary placement of the LRC alt contigs.

    Intervals that do are the ones whose MAPQ-20 counts collapse in a
    non-ALT-aware alignment.
    """
    c0, s0, e0 = ALT_PLACEMENT
    return chrom == c0 and start < e0 and end > s0


def gene_span(name: str) -> tuple[str, int, int]:
    """Primary-assembly interval of a gene body.

    Raises for LILRA3, which has none — callers that could reach it with LILRA3
    should be using :data:`LILRA3_ALT` and know why.
    """
    g = GENE_BODIES[name]
    if not g.on_primary:
        raise ValueError(
            f"{name} is not on the {ASSEMBLY} primary assembly ({g.note}); "
            "use LILRA3_ALT and the junction assay instead"
        )
    return (g.chrom, g.start, g.end)


def extraction_regions(include_alts: bool = True) -> list[str]:
    """Region strings for ``samtools view``, covering everything a LILR read
    could be aligned to: the primary cluster, and the alt contigs that carry
    LILRA3 and its neighbours.
    """
    chrom, start, end = LRC_SLICE
    regions = [f"{chrom}:{start + 1}-{end}"]   # samtools regions are 1-based
    if include_alts:
        regions.extend(ALT_CONTIGS)
    return regions
