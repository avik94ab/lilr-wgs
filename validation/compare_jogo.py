#!/usr/bin/env python3
"""Score this pipeline's LILRA6 copy number against JoGo-LILR's published calls.

The HPRC overlap scores against assemblies; this scores against another caller.
JoGo-LILR (Nagasaki et al., *Hum Immunol* 2025;86(3):111272) is cohort-relative
and reads the LILRB3+LILRA6 pair total rather than the paralogue-unique window,
so it fails differently — which is the whole value of the comparison. Its release
ships the paper's own calls for 3,202 1000 Genomes samples, so any 1KGP cohort
run here is already in it:

    JoGo-LILR_v1.1/test/paper.hapmap3202.allele_stable.tsv

    python validation/compare_jogo.py results/kgp2504/cn_calls.tsv \\
        --published $JOGO/test/paper.hapmap3202.allele_stable.tsv \\
        -o validation/reports/jogo_kgp2504.txt

**Their call is a haplotype pair, ours is a diploid count**, so the comparison
runs on what both can express. `CNV4_B1A0_B1A2` means two haplotypes carrying
B1A0 and B1A2 — diploid LILRA6 2. Where their string holds several alternative
decompositions (`CNV4_B1A0_B1A2_CNV4_B1A1_B1A1`) they are alternative phasings of
the same diploid count, and a sample whose alternatives disagree on that count is
reported as unparseable rather than guessed at.

**LILRB3 is deliberately not scored.** Their `diploid_stable` anchor grid contains
only B1/B2 haplotype types, so LILRB3 CN 2 is very nearly the only answer the
method can return; agreeing with it would be close to agreeing that most people
have two copies of LILRB3, which both methods manage without being right about
anything. See `validation/jogo_crosscheck.md`.
"""

from __future__ import annotations

import argparse
import collections
import csv
import re
import sys
from pathlib import Path

HAPLOTYPE = re.compile(r"B(\d+)A(\d+)")


def lilra6_from_type(allele_type: str) -> int | None:
    """Diploid LILRA6 copy number from a `CNV4_B1A0_B1A2`-style string.

    None when the alternative decompositions disagree, which is the one case
    where their call does not determine a diploid count.
    """
    counts = {
        sum(int(a) for _, a in HAPLOTYPE.findall(solution))
        for solution in allele_type.split("_CNV") if HAPLOTYPE.search(solution)
    }
    return counts.pop() if len(counts) == 1 else None


def read_published(path: Path) -> dict[str, int]:
    published = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            call = lilra6_from_type(row["selected_cnv_type"])
            if call is not None:
                published[row["sampleid"]] = call
    return published


def read_calls(path: Path, gene: str = "LILRA6") -> dict[str, dict]:
    with open(path) as fh:
        return {row["sample"]: row for row in csv.DictReader(fh, delimiter="\t")
                if row["gene"] == gene}


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("calls", type=Path, help="cn_calls.tsv from this pipeline")
    p.add_argument("--published", type=Path, required=True,
                   help="paper.hapmap3202.allele_stable.tsv from the JoGo release")
    p.add_argument("-o", "--output", type=Path)
    args = p.parse_args()

    published = read_published(args.published)
    ours = read_calls(args.calls)
    shared = sorted(set(ours) & set(published))

    agree = 0
    disagreements = []
    confusion: collections.Counter = collections.Counter()
    # A refusal is not a disagreement and not an agreement: it is the pipeline
    # declining to call, and collapsing it into either would hide the thing the
    # status column exists to show.
    refused = [s for s in shared if ours[s]["status"] != "measured"]
    for sample in shared:
        row = ours[sample]
        if row["status"] != "measured":
            continue
        called, expected = int(row["copies"]), published[sample]
        confusion[(expected, called)] += 1
        if called == expected:
            agree += 1
        else:
            disagreements.append((sample, expected, called, row["estimate"],
                                  row["confidence"], row["ambiguous"]))

    scored = len(shared) - len(refused)
    lines = ["LILRA6 copy number vs JoGo-LILR published calls",
             f"  our calls          {len(ours)}",
             f"  their published    {len(published)}",
             f"  in both            {len(shared)}",
             f"  refused by us      {len(refused)}"
             " (status != measured; not scored either way)",
             ""]
    if scored:
        lines.append(f"agreement {agree}/{scored} = {100 * agree / scored:.1f}%")
    lines.append("")

    lines.append("confusion (their CN down, ours across)")
    values = sorted({cn for pair in confusion for cn in pair})
    lines.append("      " + "".join(f"{v:>7}" for v in values))
    for expected in values:
        row = "".join(f"{confusion.get((expected, called), 0):>7}" for called in values)
        lines.append(f"{expected:>5} {row}")
    lines.append("")

    if disagreements:
        lines.append("disagreements")
        lines.append(f"  {'sample':<10}{'theirs':>7}{'ours':>6}{'estimate':>10}"
                     f"{'conf':>7}  flagged")
        for sample, expected, called, estimate, conf, flagged in disagreements:
            lines.append(f"  {sample:<10}{expected:>7}{called:>6}{estimate:>10}"
                         f"{conf:>7}  {flagged}")
        lines.append("")

    lines.append("Read with the caveats in validation/jogo_crosscheck.md: LILRB3 is not")
    lines.append("scored here because their anchor grid makes agreement there nearly")
    lines.append("automatic, and both methods are being compared on a cohort drawn from")
    lines.append("the same 1000 Genomes collection their background table is built from.")

    report = "\n".join(lines)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(report + "\n")
    print(report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
