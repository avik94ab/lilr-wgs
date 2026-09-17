#!/usr/bin/env python3
"""Collect every sample's coverage model into one table.

The ALT-awareness verdict is the reason this is its own output rather than a
column buried in the summary. One sample reading `not_alt_aware` is a curiosity;
the whole cohort reading it means no MAPQ-20 window in the dataset means
anything, and every LILRA6 and LILRB3 call in it is "not measured" rather than a
number. That is a fact about the dataset, so it belongs somewhere a person will
actually look.
"""

from __future__ import annotations

import csv
import glob
import json
import statistics
import sys
from collections import Counter
from pathlib import Path

COLUMNS = ["sample", "lambda1", "efficiency", "n_efficiency_anchors",
           "dispersion", "lambda1_outside", "q20_lrc", "q20_outside",
           "dilution", "alt_verdict", "usable_mapq20", "n_controls", "warnings"]


def main() -> int:
    cov_dir, out_path = sys.argv[1:3]
    rows = []
    for path in sorted(glob.glob(f"{cov_dir}/*.json")):
        try:
            rows.append(json.loads(Path(path).read_text()))
        except Exception:
            print(f"warning: could not read {path}", file=sys.stderr)

    if not rows:
        raise SystemExit(f"no coverage models under {cov_dir}")

    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    verdicts = Counter(r.get("alt_verdict", "?") for r in rows)
    lambdas = [float(r["lambda1"]) for r in rows if r.get("lambda1")]
    effs = [float(r["efficiency"]) for r in rows if r.get("efficiency")]

    print(f"wrote {out} ({len(rows)} samples)")
    if lambdas:
        print(f"  lambda1:    median {statistics.median(lambdas):.2f} "
              f"range {min(lambdas):.1f}-{max(lambdas):.1f}")
    if effs:
        print(f"  efficiency: median {statistics.median(effs):.3f} "
              f"range {min(effs):.3f}-{max(effs):.3f}")
    print(f"  verdicts:   {dict(verdicts)}")

    bad = sum(n for v, n in verdicts.items() if v != "alt_aware")
    if bad:
        print(f"  WARNING: {bad}/{len(rows)} samples are not cleanly ALT-aware. "
              "LILRA6 and LILRB3 are reported as not measured in those, rather "
              "than as zero.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
