#!/usr/bin/env python3
"""Score realigned LILRA6 copy number against the calls made from the CRAM as-is.

    python3 validation/compare_realign.py \
        --realigned results/lilra6/lilra6_cn.tsv \
        --as-is results/kgp100/cn/cn_calls.tsv

Both sides measure the same reads from the same samples. What differs is one
step: the as-is calls read depth out of the alignment NYGC produced, and the
realigned ones pull those reads back to FASTQ and align them again, against the
same index, with this pipeline's own parameters.

**Agreement here is a plumbing test, not a validation.** Because the realignment
targets the index the input was aligned against, near-identical placement is the
expected result, and the thing it rules out is that extraction, collation,
re-pairing, singleton handling or the `-Y`/`-K` settings lost or moved reads.
What it cannot show is that the pipeline now works on an input aligned to
something else -- that needs an input aligned to something else. Disagreement is
the informative direction: it localises which of those steps is lossy.

Against truth, `validation/compare_cn.py` remains the scorer; this compares two
measurements of the same thing.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


def read_calls(path: Path, gene: str) -> dict[str, dict]:
    """``{sample: row}`` for one gene, from either output's TSV layout."""
    with path.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return {r["sample"]: r for r in rows if r.get("gene") == gene}


def to_int(value: str) -> int | None:
    """An empty copy number is not zero. It is the absence of a measurement, and
    the whole point of the status column is that the two are different."""
    value = (value or "").strip()
    return int(value) if value not in ("", "NA") else None


def to_float(value: str) -> float | None:
    value = (value or "").strip()
    try:
        return float(value)
    except ValueError:
        return None


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--realigned", type=Path, required=True)
    p.add_argument("--as-is", type=Path, required=True)
    p.add_argument("--gene", default="LILRA6")
    p.add_argument("--disagreements", type=Path,
                   help="write the disagreeing samples here as a TSV")
    args = p.parse_args()

    new = read_calls(args.realigned, args.gene)
    old = read_calls(args.as_is, args.gene)
    shared = sorted(set(new) & set(old))

    print(f"{args.gene}: {len(new)} realigned, {len(old)} as-is, "
          f"{len(shared)} in common")
    only_new, only_old = sorted(set(new) - set(old)), sorted(set(old) - set(new))
    if only_new:
        print(f"  only realigned: {len(only_new)} ({', '.join(only_new[:5])} ...)")
    if only_old:
        print(f"  only as-is:     {len(only_old)} ({', '.join(only_old[:5])} ...)")
    if not shared:
        raise SystemExit("no samples in common")

    agree = 0
    both_measured = 0
    status_changed: Counter = Counter()
    deltas: list[float] = []
    disagreements: list[dict] = []

    for sample in shared:
        n, o = new[sample], old[sample]
        ns, os_ = n.get("status", ""), o.get("status", "")
        if ns != os_:
            status_changed[f"{os_} -> {ns}"] += 1

        nc, oc = to_int(n.get("copies", "")), to_int(o.get("copies", ""))
        ne, oe = to_float(n.get("estimate", "")), to_float(o.get("estimate", ""))
        if ne is not None and oe is not None:
            deltas.append(ne - oe)
        if ns == "measured" and os_ == "measured":
            both_measured += 1
            if nc == oc:
                agree += 1
            else:
                disagreements.append({
                    "sample": sample,
                    "as_is_copies": oc, "realigned_copies": nc,
                    "as_is_estimate": oe, "realigned_estimate": ne,
                    "delta": None if (ne is None or oe is None) else round(ne - oe, 3),
                    "as_is_lambda1": o.get("support", o.get("lambda1", "")),
                    "realigned_lambda1": n.get("lambda1", ""),
                })

    print(f"\nboth measured: {both_measured}/{len(shared)}")
    if both_measured:
        print(f"  integer agreement: {agree}/{both_measured} "
              f"({100 * agree / both_measured:.1f}%)")
    if status_changed:
        print("  status changes:")
        for change, count in status_changed.most_common():
            print(f"    {change}: {count}")
    else:
        print("  status changes: none")

    if deltas:
        deltas.sort()
        mid = deltas[len(deltas) // 2]
        mean = sum(deltas) / len(deltas)
        # The median shift is the number that matters: a systematic offset means
        # the realignment recovers or loses depth uniformly, which moves every
        # call toward a boundary together. Scatter around it is per-sample noise.
        print(f"\ncontinuous estimate (realigned - as-is), n={len(deltas)}:")
        print(f"  median {mid:+.4f}   mean {mean:+.4f}   "
              f"min {deltas[0]:+.4f}   max {deltas[-1]:+.4f}")
        within = sum(1 for d in deltas if abs(d) <= 0.10)
        print(f"  |delta| <= 0.10: {within}/{len(deltas)} "
              f"({100 * within / len(deltas):.1f}%)")

    if disagreements:
        print(f"\n{len(disagreements)} integer disagreements:")
        for d in disagreements[:15]:
            print(f"  {d['sample']}: as-is {d['as_is_copies']} "
                  f"({d['as_is_estimate']}) -> realigned {d['realigned_copies']} "
                  f"({d['realigned_estimate']})")
        if len(disagreements) > 15:
            print(f"  ... and {len(disagreements) - 15} more")

    if args.disagreements and disagreements:
        with args.disagreements.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(disagreements[0]),
                               delimiter="\t")
            w.writeheader()
            w.writerows(disagreements)
        print(f"\nwrote {args.disagreements}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
