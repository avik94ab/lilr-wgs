#!/usr/bin/env python3
"""Merge per-(sample, locus) genotype markers into one summary table.

One row per sample/locus/haplotype. The columns that are new relative to the
predecessor's summary are the callability ones: `callable_fraction` and the
per-reason counts. They answer the question an N-heavy consensus leaves open —
how much of the gene the data could speak to, and what stopped it where it
could not.
"""

from __future__ import annotations

import csv
import glob
import json
import sys
from collections import Counter
from pathlib import Path

COLUMNS = [
    "sample", "locus", "cn", "status", "callable_fraction",
    "n_low_depth", "n_high_depth", "n_low_mapq", "n_paralog_ambiguous",
    "phased_ok", "n_het", "n_phased", "phasing_rate",
    "haplotype", "gdna_len", "gdna_cds_len", "cdna_len", "protein_len",
    "gdna_n_pct", "protein_x_pct", "error",
]


def main() -> int:
    marker_dir, out_csv = sys.argv[1:3]
    rows, tally = [], Counter()

    for path in sorted(glob.glob(f"{marker_dir}/*.json")):
        try:
            r = json.loads(Path(path).read_text())
        except Exception:
            tally["unreadable"] += 1
            continue
        tally[r.get("status", "?")] += 1

        rate = r.get("phasing_rate")
        base = {
            "sample": r.get("sample", ""), "locus": r.get("locus", ""),
            "cn": r.get("cn", ""), "status": r.get("status", ""),
            "callable_fraction": r.get("callable_fraction", ""),
            "n_low_depth": r.get("cb_n_low_depth", ""),
            "n_high_depth": r.get("cb_n_high_depth", ""),
            "n_low_mapq": r.get("cb_n_low_mapq", ""),
            "n_paralog_ambiguous": r.get("cb_n_paralog_ambiguous", ""),
            "phased_ok": r.get("phased_ok", ""),
            "n_het": r.get("n_het", ""), "n_phased": r.get("n_phased", ""),
            "phasing_rate": f"{rate:.3f}" if isinstance(rate, (int, float)) else "",
        }
        haps = r.get("haplotypes") or []
        if haps:
            for h in haps:
                rows.append({**base, "haplotype": h["hap"],
                             "gdna_len": h["gdna_len"],
                             "gdna_cds_len": h["gdna_cds_len"],
                             "cdna_len": h["cdna_len"],
                             "protein_len": h["protein_len"],
                             "gdna_n_pct": h["gdna_n_pct"],
                             "protein_x_pct": h["protein_x_pct"], "error": ""})
        else:
            # A locus with no haplotypes still gets a row. `absent_cn0` at
            # LILRA3 is a result, not a gap, and dropping it would make a
            # deletion homozygote indistinguishable from a sample that was
            # never processed.
            rows.append({**base, "error": r.get("error", "")})

    out = Path(out_csv)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    print(f"wrote {out} ({len(rows)} rows from {sum(tally.values())} markers)")
    print("status tally:", dict(sorted(tally.items())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
