"""Assigning read pairs to genes across the paralogues.

Ported from `lilr-genotyper`'s ``filter_crossmapped.py``, whose logic is sound
and stays: each pair is aligned against all 11 HPRC pangenome gene panels, the
paired AS scores are compared, and the pair goes to the single best gene when one
wins by more than a buffer.

The part worth restating is what happens on a tie. Dropping an ambiguous pair
from every reference is only meaningful where the references are *separable*.
LILRA6 and LILRB3 share a ~4.6 kb block encoding all four Ig domains that
contains no gene-diagnostic 31-mers, and the LILRB1 cytoplasmic tail is ~92%
identical to LILRB4 — so in those blocks the tie is a property of the genes, not
of the read, and discarding on it deletes the block from both genes at once.
Measured on a 560-sample capture cohort, that cost 26% of the LILRB3 CDS, 13% of
LILRA6 and 7% of LILRB1, while the eight separable genes lost nothing. PING makes
the same choice for its one inseparable pair, KIR2DL5A/B: the competitive filter
is commented out and the two genes share a reference.

So a tie confined to one ``shared_group`` is kept in every tied gene and tagged.

**What is new here.** In the predecessor that tag died immediately: the filtered
BAMs were converted to FASTQ with ``bedtools bamtofastq``, which carries no tags,
and everything downstream treated shared-block support and gene-unique support as
identical. That is exactly the information :func:`lilrwgs.depth_model.call_position`
needs — without it the two-sided depth gate has no way to know that a shared block
legitimately carries two genes' worth of reads, and would reject ~90% of the
LILRA6 coding sequence in every sample.

The recruit-then-call structure still needs the FASTQ hop, because reads are
recruited against a 465-sequence pangenome panel and called against a single
per-locus coordinate system. So the verdict is written alongside the FASTQs as a
read-name table instead, and re-applied after realignment. The information
crosses the gap that the tag could not.
"""

from __future__ import annotations

import csv
import os
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

DEFAULT_BUFFER = 2

# Pairs that short reads cannot separate. Defaults match the predecessor's, whose
# grounds are measured rather than assumed — see the module docstring.
DEFAULT_SHARED_GROUPS = [("LILRA6", "LILRB3"), ("LILRB1", "LILRB4")]

SHARED_TAG = "ZS"           # value: how many genes the pair was tied across
SHARED_TABLE = "shared_pairs.tsv"


@dataclass
class AssignmentStats:
    """Where every pair went, which is the audit the predecessor removed.

    Its ``process_sample.sh`` piped the arbitration log through
    ``grep -vF "Ambiguous read pair"``, which hid the fact that 70.9% of LILRA6
    pairs and 62.8% of LILRB3 pairs were being discarded. Counting them is how
    that was eventually found; the counts stay.
    """

    n_unique: int = 0
    n_arbitrated: int = 0
    n_shared: int = 0
    n_dropped: int = 0
    per_gene: dict[str, int] = field(default_factory=lambda: defaultdict(int))
    per_gene_shared: dict[str, int] = field(default_factory=lambda: defaultdict(int))

    @property
    def n_total(self) -> int:
        return self.n_unique + self.n_arbitrated + self.n_shared + self.n_dropped

    def as_row(self) -> dict:
        total = max(self.n_total, 1)
        return {
            "n_pairs": self.n_total,
            "n_unique": self.n_unique,
            "n_arbitrated": self.n_arbitrated,
            "n_shared": self.n_shared,
            "n_dropped": self.n_dropped,
            "dropped_fraction": round(self.n_dropped / total, 4),
            "shared_fraction": round(self.n_shared / total, 4),
            **{f"n_{g}": n for g, n in sorted(self.per_gene.items())},
        }


def group_index(gene: str, groups: list[tuple[str, ...]]) -> int | None:
    """Which shared group a gene belongs to, or None if it is separable."""
    for i, members in enumerate(groups):
        if gene in members:
            return i
    return None


def arbitrate(
    scores: dict[str, dict[str, int]],
    *,
    buffer: int = DEFAULT_BUFFER,
    shared_groups: list[tuple[str, ...]] | None = None,
) -> tuple[dict[str, dict[str, int]], AssignmentStats]:
    """Decide, for every read pair, which gene(s) it belongs to.

    Pure: takes paired AS sums and returns assignments, so the policy can be
    tested without BAMs.

    Args:
        scores: ``{gene: {read_name: paired_AS_sum}}``.
        buffer: a gene within this much of the best score counts as tied. 2 is
            the predecessor's default; its size is largely irrelevant for the
            LILRA6/LILRB3 pair, where 99.87% of ties have a margin of exactly
            zero — the drop-on-tie *policy* was the defect, not the buffer.

    Returns:
        ``({gene: {read_name: n_tied}}, stats)``. ``n_tied`` is 1 for a pair
        assigned outright and >1 for one kept across a shared group, which is
        what reaches the depth model.
    """
    groups = shared_groups if shared_groups is not None else DEFAULT_SHARED_GROUPS
    by_read: dict[str, list[tuple[str, int]]] = defaultdict(list)
    for gene, reads in scores.items():
        for name, score in reads.items():
            by_read[name].append((gene, score))

    assigned: dict[str, dict[str, int]] = {gene: {} for gene in scores}
    stats = AssignmentStats()

    for name, hits in by_read.items():
        if len(hits) == 1:
            gene = hits[0][0]
            assigned[gene][name] = 1
            stats.n_unique += 1
            stats.per_gene[gene] += 1
            continue

        best = max(score for _, score in hits)
        passing = [(gene, score) for gene, score in hits if score >= best - buffer]

        if len(passing) == 1:
            gene = passing[0][0]
            assigned[gene][name] = 1
            stats.n_arbitrated += 1
            stats.per_gene[gene] += 1
            continue

        # Ambiguous. If every tied gene sits in one shared group, the tie is a
        # property of the genes rather than of the read: keep the pair in all of
        # them so the shared block retains coverage. A tie spanning different
        # groups is genuinely ambiguous and is dropped.
        tied_groups = {group_index(gene, groups) for gene, _ in passing}
        if len(tied_groups) == 1 and None not in tied_groups:
            for gene, _ in passing:
                assigned[gene][name] = len(passing)
                stats.per_gene[gene] += 1
                stats.per_gene_shared[gene] += 1
            stats.n_shared += 1
        else:
            stats.n_dropped += 1

    return assigned, stats


def paired_scores(bam_path: str | os.PathLike) -> dict[str, int]:
    """Paired AS sums for one gene's BAM, keyed by read name.

    Only pairs with both mates mapped are kept — a half-mapped pair carries no
    comparative information, since the score it would contribute is missing on
    one side for every gene at once.
    """
    import pysam

    mates: dict[str, dict[str, object]] = defaultdict(dict)
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue
            key = "mate1" if read.is_read1 else "mate2"
            slot = mates[read.query_name]
            if key not in slot:
                slot[key] = read

    out: dict[str, int] = {}
    for name, slot in mates.items():
        if "mate1" not in slot or "mate2" not in slot:
            continue
        try:
            out[name] = slot["mate1"].get_tag("AS") + slot["mate2"].get_tag("AS")
        except KeyError:
            out[name] = 0
    return out


def write_filtered_bams(
    bams: dict[str, str],
    assigned: dict[str, dict[str, int]],
    out_dir: str | os.PathLike,
) -> dict[str, str]:
    """Write one filtered BAM per gene, tagging shared pairs in place."""
    import pysam

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written: dict[str, str] = {}

    for gene, bam_path in bams.items():
        keep = assigned.get(gene, {})
        out_path = out_dir / f"{gene}.filtered.bam"
        with pysam.AlignmentFile(bam_path, "rb") as src:
            with pysam.AlignmentFile(str(out_path), "wb", header=src.header) as dst:
                for read in src.fetch(until_eof=True):
                    if read.is_secondary or read.is_supplementary or read.is_unmapped:
                        continue
                    n_tied = keep.get(read.query_name)
                    if n_tied is None:
                        continue
                    if n_tied > 1:
                        read.set_tag(SHARED_TAG, n_tied, value_type="i")
                    dst.write(read)
        written[gene] = str(out_path)
    return written


def write_shared_table(assigned: dict[str, dict[str, int]],
                       path: str | os.PathLike) -> int:
    """Record which pairs are shared, so the fact survives the FASTQ hop.

    Reads are recruited against a 465-sequence pangenome panel and called
    against a single per-locus reference, so there is a FASTQ in between and BAM
    tags do not cross it. This table does. Without it the depth model cannot
    distinguish a shared block — which legitimately carries two genes' worth of
    reads — from an unseparated pile-up, and would reject the former along with
    the latter.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene", "read_name", "n_tied"])
        for gene, reads in sorted(assigned.items()):
            for name, n_tied in reads.items():
                if n_tied > 1:
                    w.writerow([gene, name, n_tied])
                    n += 1
    return n


def read_shared_names(path: str | os.PathLike, gene: str) -> dict[str, int]:
    """Read names shared for one gene, as written by :func:`write_shared_table`."""
    path = Path(path)
    if not path.exists():
        return {}
    out: dict[str, int] = {}
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["gene"] == gene:
                out[row["read_name"]] = int(row["n_tied"])
    return out
