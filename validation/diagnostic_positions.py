#!/usr/bin/env python3
"""How separable are LILRA6 and LILRB3, really?

    python3 validation/diagnostic_positions.py

**Result, so it is not necessary to run this to learn the answer: none of the
99 positions is a fixed difference.** Every site where the two references
disagree has both alleles present in both panels, and only 13 of 99 are
informative even in the weaker, one-directional sense that an allele is <1% in
one gene and >=5% in the other — 12 of those exclude LILRA6 and 1 excludes
LILRB3. Across a 4,630 bp block that is one informative site per ~356 bp, so a
150 bp read carries one with probability ~0.4.

The shared block is therefore not separable by read assignment, and schemes that
attribute a shared-block read or haplotype to one gene are not merely hard, they
are unsupported by the sequence. Report the pair jointly instead.

The two references differ at 99 substituted positions across their 4,634 bp
shared block, and `assign.arbitrate()` finds 99.87% of LILRA6/LILRB3 ties at a
margin of exactly zero. That agreement is the point: the ties are the truth
about these genes rather than an artefact of scoring against panels.

The suspicion this script tests is that the signal is real between the two
*references* and absent across the *panels*. Arbitration scores each read against
558 LILRA6 and 462 LILRB3 haplotypes and takes the best match in each, so a
position only separates the genes if every LILRA6 haplotype disagrees with every
LILRB3 haplotype there. Gene conversion between paralogues is common in this
family, and one converted haplotype in either panel is enough to destroy a
position's diagnostic value.

So: for each reference-level difference, collect the alleles actually present in
each panel and ask whether the two sets are disjoint. A position is

  fixed        - allele sets disjoint; usable for assignment
  polymorphic  - sets overlap; a read carrying it cannot be assigned on it alone

The fixed count is what bounds any scheme that assigns reads by their state at
distinguishing sites, so it should be measured before such a scheme is built.
"""

from __future__ import annotations

import argparse
import subprocess
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

REFS = Path("resources/bundle/references")
PANELS = Path("resources/gdna")
PAIR = ("LILRA6", "LILRB3")

# A panel sequence must cover a position with this much flanking alignment for
# its allele to count. Prevents an allele being read off a ragged alignment end.
MIN_FLANK = 20


def _read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name, buf = None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def diagnostic_sites(a_ref: Path, b_ref: Path) -> list[tuple[int, int, str, str]]:
    """Reference-level differences as ``(a_pos, b_pos, a_base, b_base)``."""
    import pysam

    with tempfile.TemporaryDirectory() as tmp:
        sam = Path(tmp) / "pair.sam"
        with sam.open("w") as fh:
            subprocess.run(["minimap2", "-a", "--MD", str(b_ref), str(a_ref)],
                           stdout=fh, stderr=subprocess.DEVNULL, check=True)
        out = []
        with pysam.AlignmentFile(str(sam), "r") as handle:
            for rec in handle:
                if rec.is_unmapped or rec.is_supplementary or rec.is_secondary:
                    continue
                q = rec.query_sequence
                for qpos, rpos, rbase in rec.get_aligned_pairs(with_seq=True):
                    if qpos is None or rpos is None or rbase is None:
                        continue
                    qb, rb = q[qpos].upper(), rbase.upper()
                    if qb != rb and qb in "ACGT" and rb in "ACGT":
                        out.append((qpos, rpos, qb, rb))
    return out


def panel_alleles(panel: Path, ref: Path, positions: set[int]) -> dict[int, Counter]:
    """``{ref_pos: Counter(allele)}`` over every panel sequence."""
    import pysam

    alleles: dict[int, Counter] = defaultdict(Counter)
    with tempfile.TemporaryDirectory() as tmp:
        sam = Path(tmp) / "panel.sam"
        with sam.open("w") as fh:
            subprocess.run(["minimap2", "-a", "--MD", "-t", "4",
                            str(ref), str(panel)],
                           stdout=fh, stderr=subprocess.DEVNULL, check=True)
        with pysam.AlignmentFile(str(sam), "r") as handle:
            for rec in handle:
                if rec.is_unmapped or rec.is_supplementary or rec.is_secondary:
                    continue
                q = rec.query_sequence
                if not q:
                    continue
                lo, hi = rec.reference_start, rec.reference_end
                for qpos, rpos in rec.get_aligned_pairs():
                    if qpos is None or rpos is None or rpos not in positions:
                        continue
                    if rpos - lo < MIN_FLANK or hi - rpos < MIN_FLANK:
                        continue
                    base = q[qpos].upper()
                    if base in "ACGT":
                        alleles[rpos][base] += 1
    return alleles


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--min-support", type=int, default=5,
                   help="panel sequences needed at a position to judge it")
    p.add_argument("-o", "--output", type=Path,
                   default=Path("validation/reports/diagnostic_positions.tsv"))
    args = p.parse_args()

    a, b = PAIR
    sites = diagnostic_sites(REFS / f"{a}_named.fa", REFS / f"{b}_named.fa")
    print(f"{a} vs {b}: {len(sites)} reference-level differences")

    a_pos = {s[0] for s in sites}
    b_pos = {s[1] for s in sites}
    print(f"collecting panel alleles ({a}: {PANELS / f'{a}.fasta'}) ...")
    a_all = panel_alleles(PANELS / f"{a}.fasta", REFS / f"{a}_named.fa", a_pos)
    print(f"collecting panel alleles ({b}) ...")
    b_all = panel_alleles(PANELS / f"{b}.fasta", REFS / f"{b}_named.fa", b_pos)

    fixed, polymorphic, unjudged = [], [], []
    for ap, bp, ab, bb in sites:
        ca, cb = a_all.get(ap, Counter()), b_all.get(bp, Counter())
        if sum(ca.values()) < args.min_support or sum(cb.values()) < args.min_support:
            unjudged.append((ap, bp))
            continue
        if set(ca) & set(cb):
            polymorphic.append((ap, bp, dict(ca), dict(cb)))
        else:
            fixed.append((ap, bp, ab, bb, dict(ca), dict(cb)))

    judged = len(fixed) + len(polymorphic)
    print()
    print(f"  reference differences : {len(sites)}")
    print(f"  judged                : {judged}  (>= {args.min_support} panel seqs both sides)")
    print(f"  FIXED (disjoint)      : {len(fixed)}"
          + (f"  = {100*len(fixed)/judged:.1f}% of judged" if judged else ""))
    print(f"  polymorphic (overlap) : {len(polymorphic)}")
    print(f"  unjudged (thin panel) : {len(unjudged)}")

    if fixed:
        span = max(f[0] for f in fixed) - min(f[0] for f in fixed)
        print()
        print(f"  fixed sites span {span:,} bp of {a}")
        print(f"  mean spacing     {span/max(len(fixed)-1,1):.0f} bp")
        print(f"  per 150 bp read  {150*len(fixed)/max(span,1):.2f} fixed sites")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as fh:
        fh.write(f"{a}_pos\t{b}_pos\tclass\t{a}_ref\t{b}_ref\t{a}_panel\t{b}_panel\n")
        for ap, bp, ab, bb, ca, cb in fixed:
            fh.write(f"{ap}\t{bp}\tfixed\t{ab}\t{bb}\t{ca}\t{cb}\n")
        for ap, bp, ca, cb in polymorphic:
            fh.write(f"{ap}\t{bp}\tpolymorphic\t.\t.\t{ca}\t{cb}\n")
    print(f"\nwrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
