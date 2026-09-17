#!/usr/bin/env python3
"""Score called copy number against the HPRC truth set.

Reports exact-match accuracy and the confusion matrix per gene, and — because
this is the number that decides whether the pipeline is usable — reports them
separately for calls the pipeline flagged as ambiguous and calls it did not. A
method that is 95% accurate overall but 99% accurate on its confident calls and
60% on its flagged ones is a different, and much more useful, method than one
that is uniformly 95%.

Two caveats are printed rather than buried, because forgetting either would
overstate the result:

* **Circularity.** A donor in the overlap is aligned against a panel containing
  its own haplotypes. Pass `--leave-one-donor-out` results for the number that
  generalises.
* **Inferred absences.** LILRA3 CN 0 in the truth set comes from a donor being
  absent from that panel while present in the others, not from counting. Use
  `--counted-only` to score without them.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


def load_truth(path: Path, counted_only: bool) -> dict[tuple[str, str], int]:
    out = {}
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if counted_only and row.get("source") == "inferred_absent":
                continue
            out[(row["donor"], row["gene"])] = int(row["copies"])
    return out


def load_calls(path: Path) -> dict[tuple[str, str], dict]:
    out = {}
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[(row["sample"], row["gene"])] = row
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("calls", type=Path, help="results/cn/cn_calls.tsv")
    p.add_argument("--truth", type=Path,
                   default=Path("validation/truth/copy_number.tsv"))
    p.add_argument("--counted-only", action="store_true",
                   help="exclude truth entries inferred from panel absence")
    p.add_argument("--leave-one-donor-out", action="store_true",
                   help="label the report as the generalising run")
    p.add_argument("-o", "--output", type=Path)
    args = p.parse_args()

    truth = load_truth(args.truth, args.counted_only)
    calls = load_calls(args.calls)

    per_gene: dict[str, Counter] = defaultdict(Counter)
    confusion: dict[str, Counter] = defaultdict(Counter)
    for key, row in calls.items():
        if key not in truth:
            continue
        gene = key[1]
        status = row.get("status", "")
        if status != "measured":
            per_gene[gene][f"not_scored_{status}"] += 1
            continue
        called = int(row["copies"])
        expected = truth[key]
        flagged = row.get("ambiguous", "").lower() in ("true", "1")
        bucket = "flagged" if flagged else "confident"
        per_gene[gene][f"{bucket}_n"] += 1
        if called == expected:
            per_gene[gene][f"{bucket}_correct"] += 1
        confusion[gene][(expected, called)] += 1

    lines = []
    mode = ("leave-one-donor-out" if args.leave_one_donor_out
            else "AS-IS (donor's own haplotypes are in the panel; "
                 "this overstates accuracy)")
    lines.append(f"copy-number accuracy vs HPRC truth — {mode}")
    if not args.counted_only:
        lines.append("truth includes CN 0 inferred from panel absence "
                     "(--counted-only to exclude)")
    lines.append("")
    header = f"{'gene':8s} {'n':>5s} {'confident':>12s} {'flagged':>12s} {'overall':>9s}"
    lines.append(header)

    for gene in sorted(per_gene):
        c = per_gene[gene]
        cn_, cc = c["confident_n"], c["confident_correct"]
        fn, fc = c["flagged_n"], c["flagged_correct"]
        total, correct = cn_ + fn, cc + fc
        if total == 0:
            continue
        lines.append(
            f"{gene:8s} {total:5d} "
            f"{(f'{cc}/{cn_} {cc / cn_:.1%}' if cn_ else '-'):>12s} "
            f"{(f'{fc}/{fn} {fc / fn:.1%}' if fn else '-'):>12s} "
            f"{correct / total:9.1%}")

    lines.append("")
    lines.append("confusion (truth -> called), where they differ:")
    for gene in sorted(confusion):
        wrong = {k: v for k, v in confusion[gene].items() if k[0] != k[1]}
        if wrong:
            detail = ", ".join(f"{t}->{c}: {n}" for (t, c), n in sorted(wrong.items()))
            lines.append(f"  {gene}: {detail}")

    skipped = {g: {k: v for k, v in c.items() if k.startswith("not_scored")}
               for g, c in per_gene.items()}
    skipped = {g: v for g, v in skipped.items() if v}
    if skipped:
        lines.append("")
        lines.append("not scored (the pipeline declined to call):")
        for gene, counts in sorted(skipped.items()):
            lines.append(f"  {gene}: {dict(counts)}")

    report = "\n".join(lines)
    print(report)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(report + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
