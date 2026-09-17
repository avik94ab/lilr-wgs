#!/usr/bin/env python3
"""Fetch every sample's LRC slice to local disk, and emit a manifest for them.

For sites where the machine that can reach the internet is not the machine that
runs the jobs. On Wynton the compute nodes have no outbound route at all — a
remote CRAM there fails with "Destination address required" — so the pipeline's
one remote pass has to happen somewhere else and land on shared storage first.

    python3 scripts/stage_slices.py --manifest config/manifest.hprc.tsv \
        --reference resources/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa \
        --outdir /wynton/scratch/$USER/lilr_slices \
        --out-manifest config/manifest.hprc.staged.tsv --jobs 6

The slice is the designed staging point rather than a convenient one: it already
carries every region any stage needs — the LRC, the LILRA3 alt intervals, the
deletion junction, and the control loci both inside and outside the alt
placement — because `loci.slice_intervals` collects them for exactly this reason.
So a staged slice is a drop-in for the CRAM it came from, and `process_sample`
takes it with no flag: re-slicing an already-sliced BAM is a local no-op.

Costs ~12 MB and ~20 s of network per sample, and it is resumable — a sample
whose slice is already present and non-empty is skipped, so an interrupted run
costs only what it had not finished.

Deliberately network-bound and nothing else: no coverage model, no copy number.
Those are CPU, they belong on the cluster, and keeping them out of here means the
staging step can be run on a login node without being antisocial.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from lilrwgs import extract  # noqa: E402
from lilrwgs.shell import ToolError  # noqa: E402


def stage_one(row: dict, reference: Path, outdir: Path, threads: int,
              tmpdir: str | None, force: bool) -> dict:
    """One sample's slice. Returns a status row; never raises."""
    sample = row["sample_id"]
    out_bam = outdir / f"{sample}.slice.bam"

    if out_bam.exists() and out_bam.stat().st_size > 0 and not force:
        return {"sample": sample, "bam": str(out_bam), "status": "cached"}

    try:
        info = extract.slice_cram(sample, row["cram"], reference, out_bam,
                                  threads=threads, tmpdir=tmpdir)
        return {**info, "status": "staged"}
    except (ToolError, OSError) as exc:
        # A failed fetch must not leave a truncated BAM behind: the next run
        # would treat it as cached, and a short slice is a low copy number
        # rather than an error.
        for leftover in (out_bam, out_bam.with_suffix(".bam.bai")):
            leftover.unlink(missing_ok=True)
        return {"sample": sample, "status": "failed", "error": str(exc)[:400]}


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--manifest", required=True,
                   help="TSV with sample_id, cram, crai, population")
    p.add_argument("--reference", required=True)
    p.add_argument("--outdir", required=True, help="where the slices land")
    p.add_argument("--out-manifest", required=True,
                   help="manifest rewritten to point at the staged slices")
    p.add_argument("--jobs", type=int, default=6,
                   help="concurrent fetches. This is network-bound, so more than "
                        "the core count is fine and more than ~8 mostly annoys EBI")
    p.add_argument("--threads", type=int, default=2, help="samtools threads per fetch")
    p.add_argument("--tmpdir", default=os.environ.get("TMPDIR"))
    p.add_argument("--force", action="store_true", help="restage even if present")
    p.add_argument("--stats", help="write per-sample staging stats here as JSON")
    args = p.parse_args(argv)

    reference = Path(args.reference).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    with open(args.manifest) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        print(f"{args.manifest} has no samples", file=sys.stderr)
        return 1

    extract.check_remote_support()

    started = time.time()
    results: dict[str, dict] = {}
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {
            pool.submit(stage_one, r, reference, outdir, args.threads,
                        args.tmpdir, args.force): r["sample_id"]
            for r in rows
        }
        for i, fut in enumerate(as_completed(futures), 1):
            res = fut.result()
            results[res["sample"]] = res
            note = res.get("error", f"{res.get('n_records', '-')} records")
            print(f"[{i}/{len(rows)}] {res['sample']}: {res['status']} ({note})",
                  flush=True)

    ok = [r for r in rows if results[r["sample_id"]]["status"] != "failed"]
    failed = [s for s, r in results.items() if r["status"] == "failed"]

    # Only samples that actually have a slice go into the manifest. A row
    # pointing at a missing file would fail later, further from the cause.
    out_manifest = Path(args.out_manifest)
    out_manifest.parent.mkdir(parents=True, exist_ok=True)
    with out_manifest.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample_id", "cram", "crai", "population"],
                           delimiter="\t")
        w.writeheader()
        for r in ok:
            bam = outdir / f"{r['sample_id']}.slice.bam"
            w.writerow({"sample_id": r["sample_id"], "cram": str(bam),
                        "crai": str(bam) + ".bai",
                        "population": r.get("population", "")})

    no_alts = [s for s, r in results.items() if r.get("has_alt_contigs") is False]
    if no_alts:
        print(f"warning: {len(no_alts)} CRAMs have no LRC alt contigs; LILRA3 "
              f"is not measurable by depth in those: {', '.join(no_alts[:5])}",
              file=sys.stderr)

    if args.stats:
        Path(args.stats).write_text(json.dumps(results, indent=2, sort_keys=True))

    print(f"\nstaged {len(ok)}/{len(rows)} samples in "
          f"{round(time.time() - started)}s -> {out_manifest}")
    if failed:
        print(f"failed: {', '.join(failed)}", file=sys.stderr)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
