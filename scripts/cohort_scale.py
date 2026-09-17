#!/usr/bin/env python3
"""Merge per-sample copy-number calls and check the cohort's scale.

A leaf in the DAG, not a barrier. With an absolute baseline the per-copy unit
should already be 1.0; if the cohort's estimates cluster on something else,
lambda_1 is systematically off by that factor and that is a finding to
investigate, not something to divide out. Nothing downstream consumes this.
"""

from __future__ import annotations

import csv
import glob
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from lilrwgs.cn import CNCall, refine_cohort  # noqa: E402


def main() -> int:
    cn_dir, merged_path, report_path = sys.argv[1:4]

    rows, calls = [], []
    for path in sorted(glob.glob(f"{cn_dir}/*.tsv")):
        with open(path) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                rows.append(row)
                calls.append(CNCall(
                    sample=row["sample"], gene=row["gene"],
                    estimate=float(row["estimate"]) if row["estimate"] else None,
                    copies=int(row["copies"]) if row["copies"] not in ("", None) else None,
                    confidence=float(row["confidence"] or 0),
                    method=row.get("method", ""), status=row.get("status", ""),
                ))

    if not rows:
        raise SystemExit(f"no per-sample CN files under {cn_dir}")

    merged = Path(merged_path)
    merged.parent.mkdir(parents=True, exist_ok=True)
    with merged.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    report = refine_cohort(calls)
    # Copy-number distribution per gene, restricted to measured calls. Reported
    # per status so that a gate firing unevenly across copy-number classes is
    # visible: such a gate biases the allele frequency, and frequencies should
    # be computed on the unfiltered calls.
    for gene in sorted({c.gene for c in calls}):
        gene_calls = [c for c in calls if c.gene == gene]
        report.setdefault(gene, {})
        report[gene]["status_counts"] = dict(Counter(c.status for c in gene_calls))
        report[gene]["cn_distribution"] = dict(sorted(Counter(
            c.copies for c in gene_calls if c.status == "measured").items()))

    Path(report_path).write_text(json.dumps(report, indent=2))
    print(f"wrote {merged} ({len(rows)} calls) and {report_path}")
    for gene, info in sorted(report.items()):
        print(f"  {gene}: {info.get('note', '')} "
              f"cn={info.get('cn_distribution', {})}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
