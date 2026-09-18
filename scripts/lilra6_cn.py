#!/usr/bin/env python3
"""LILRA6 copy number from a list of CRAMs, by extracting and realigning.

    python3 scripts/lilra6_cn.py --inputs crams.tsv \
        --reference resources/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa \
        --threads 8 --jobs 4 -o lilra6_cn.tsv

`--inputs` is a list file, one sample per line, whitespace- or tab-separated:

    /data/HG00096.final.cram    /data/HG00096.final.cram.crai
    https://.../HG00097.final.cram
    HG00099  /data/odd/name.cram  /indexes/odd.crai

One, two or three columns. With one, the index is wherever htslib finds it by
suffix and the sample name comes from the filename; with two, the second column
is the index; with three, the first is an explicit sample name. A `#` line is a
comment, and a header line naming the columns is tolerated and skipped.

Why this exists rather than `process_sample.py`: that one reads the depth
straight out of the input alignment, so every copy number it reports is
conditional on the input having been aligned to GRCh38 *with the `.alt` file*.
That holds for the 1000 Genomes 30x set and is an assumption everywhere else.
Here the LRC and the control loci are pulled out of whatever alignment they
arrived in and realigned against the analysis set with its alt index, so
ALT-awareness is a property of this pipeline rather than of the input.

What comes out is LILRA6 and only LILRA6, one row per sample.

`status` is not decoration: a failed measurement and a true zero are different
values. LILRA6 CN 0 is real and rare, and if the realignment is not ALT-aware
every MAPQ-20 window in the cluster reads near zero for every sample alike — so
that case is reported `not_measured` and left empty, never filled in as 0.

LILRB3 is still measured internally, because LILRA6's only independent check is
the pooled LILRA6+LILRB3 depth; it is not reported. LILRA3 is not measured at
all. Its depth route cannot survive a regional extraction — the reads it counts
have their primaries scattered genome-wide — and its junction route, while
usable, is materially weaker here than on a CRAM slice. `cn.call_sample` is
where LILRA3 is called, from the CRAM as-is. See PLAN.md §12.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from lilrwgs import cn, coverage, realign  # noqa: E402
from lilrwgs.shell import ToolError  # noqa: E402

# Stripped from a filename to get a sample name, longest first so that
# `.final.cram` wins over `.cram` and does not leave a trailing `.final`.
SAMPLE_SUFFIXES = (".final.cram", ".slice.bam", ".cram", ".bam", ".sam")

# Column names a header line might use, so one can be recognised and skipped
# rather than silently processed as a sample called "sample_id".
HEADER_TOKENS = {"sample", "sample_id", "cram", "crai", "bam", "bai", "index"}

OUTPUT_FIELDS = [
    "sample", "gene", "copies", "estimate", "confidence", "status",
    "ambiguous", "method", "lambda1", "mean_depth", "alt_verdict",
    "assembly", "n_pairs", "notes",
    # The raw counts the call was made from. Carried for the same reason
    # `cn.CNCall.as_row` carries it: a copy number without its evidence cannot be
    # re-adjudicated later, and in this cluster the calls that need
    # re-adjudicating are exactly the plausible ones. Dropping this column cost
    # an afternoon -- a systematic LILRA3 shift showed up in the comparison and
    # the numbers needed to explain it had been thrown away.
    "support",
]


def sample_name(path: str) -> str:
    base = os.path.basename(path.split("?")[0].rstrip("/"))
    for suffix in SAMPLE_SUFFIXES:
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return os.path.splitext(base)[0]


def read_inputs(path: Path) -> list[dict]:
    """The list file -> ``[{sample, source, index}]``.

    Deliberately forgiving about shape and strict about duplicates: a list
    assembled by hand or by `ls` is the normal case, but two rows claiming the
    same sample name means one of them silently overwrites the other's output.
    """
    rows: list[dict] = []
    seen: dict[str, str] = {}

    for lineno, raw in enumerate(path.read_text().splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t") if "\t" in line else line.split()
        fields = [f for f in (f.strip() for f in fields) if f]
        if not fields:
            continue
        if lineno == 1 and {f.lower() for f in fields} & HEADER_TOKENS:
            continue

        if len(fields) == 1:
            sample, source, index = sample_name(fields[0]), fields[0], None
        elif len(fields) == 2:
            sample, source, index = sample_name(fields[0]), fields[0], fields[1]
        else:
            sample, source, index = fields[0], fields[1], fields[2]

        if sample in seen:
            raise SystemExit(
                f"{path}:{lineno}: two rows both name the sample {sample!r}\n"
                f"  {seen[sample]}\n  {source}\n"
                "Give the first column an explicit, distinct sample name."
            )
        seen[sample] = source
        rows.append({"sample": sample, "source": source, "index": index})

    if not rows:
        raise SystemExit(f"{path}: no inputs found")
    return rows


def genes_for(target_assembly: str, requested: str | None = None,
              ) -> tuple[str, ...]:
    """Which genes this target can report.

    Used for the failure path as well as the success one, so a sample that dies
    leaves the same set of rows it would have produced. A gene silently missing
    for some samples and present for others is a third state nobody reads for.
    """
    mod = realign.LOCI_BY_ASSEMBLY.get(target_assembly)
    windows = getattr(mod, "UNIQUE_WINDOWS", {}) if mod else {}
    available = ("LILRA6", "LILRA3") if "LILRA3" in windows else ("LILRA6",)
    if requested is None:
        return available
    asked = tuple(g.strip().upper() for g in requested.split(",") if g.strip())
    missing = [g for g in asked if g not in available]
    if missing:
        raise SystemExit(
            f"{target_assembly} cannot report {', '.join(missing)}; "
            f"it has windows for {', '.join(available)}")
    return asked


def call_one(row: dict, reference: str, bwa_index: str, outdir: Path,
             threads: int, keep_bam: bool, target_assembly: str = "GRCh38",
             target: str | None = None,
             genes: tuple[str, ...] = ("LILRA6",)) -> dict:
    """One sample, end to end. Returns a status dict; never raises.

    Never raises because a cohort is a list of independent samples and one
    unreadable CRAM should cost that sample, not the run. The failure is carried
    in the row rather than in the exit code, and `status` is `failed` — which is
    a different value from a copy number of 0 everywhere downstream.
    """
    sample = row["sample"]
    work = Path(os.environ.get("TMPDIR", "/tmp")) / f"lilra6_{sample}"
    work.mkdir(parents=True, exist_ok=True)
    bam = (outdir / "realigned" / f"{sample}.bam") if keep_bam else \
        (work / f"{sample}.bam")
    started = time.time()

    try:
        stats = realign.realign_sample(
            sample, row["source"], reference, bwa_index, bam,
            threads=threads, tmpdir=str(work), index=row["index"])

        # The coverage model reads the realigned BAM, so lambda_1 and the LILRA6
        # window depth are measured in the same units on the same alignment --
        # which is what makes the ratio a copy number rather than a comparison
        # between two pipelines' losses.
        #
        # Both the coordinate table and the FASTA are the *target's*, not the
        # input's. The table because the BAM is in the target's coordinates; the
        # FASTA because the GC correction reads sequence at the control loci to
        # build the curve, and reading it from the wrong assembly would correct
        # λ₁ by the GC of the wrong 3 Mb. Passing None there does not error --
        # it silently drops the correction, which config.yaml turns on because
        # the LILR genes are not at the genomic mean GC.
        target_loci = realign.LOCI_BY_ASSEMBLY[target_assembly]
        model = coverage.measure(sample, str(bam), reference=target or bwa_index,
                                 loci_mod=target_loci)
        calls = []
        if "LILRA6" in genes:
            calls.append(cn.call_lilra6(sample, str(bam), model,
                                        reference=target or bwa_index,
                                        loci_mod=target_loci))
        # LILRA3 only where the assembly carries it. On GRCh38 this returns
        # `not_measured` with a reason rather than a number, because the routes
        # that do work there -- alt-contig depth, the deletion junction -- need
        # the CRAM slice, not a realigned regional extraction.
        if "LILRA3" in genes:
            calls.append(cn.call_lilra3_primary(
                sample, str(bam), model, reference=target or bwa_index,
                loci_mod=target_loci))
        rows = [_row(sample, c, model, stats, target_assembly) for c in calls]
        return {"sample": sample, "ok": True, "rows": rows,
                "realign": stats.as_row(), "coverage": model.as_row(),
                "elapsed_s": round(time.time() - started, 1)}

    except (ToolError, OSError, ValueError) as exc:
        # A row, not a gap. A sample missing from the output is indistinguishable
        # from one nobody asked for; `failed` with an empty `copies` is a third
        # thing, distinct from both a measurement and a zero.
        return {
            "sample": sample, "ok": False, "error": str(exc)[:600],
            "elapsed_s": round(time.time() - started, 1),
            "rows": [{"sample": sample, "gene": g, "copies": "", "estimate": "",
                      "confidence": 0.0, "status": "failed", "ambiguous": False,
                      "method": "", "lambda1": "", "mean_depth": "",
                      "alt_verdict": "", "assembly": "", "n_pairs": "",
                      "support": "",
                      "notes": str(exc)[:200].replace("\n", " ")}
                     for g in genes],
        }
    finally:
        # Always: when --keep-realigned is set the BAM is written under outdir,
        # not here, so there is nothing in the scratch dir worth keeping either
        # way. Left behind, it is ~200 MB of FASTQ and slice per sample.
        shutil.rmtree(work, ignore_errors=True)


def _row(sample: str, call, model, stats, target_assembly: str) -> dict:
    """One output row from one CNCall."""
    r = call.as_row()
    return {
        "sample": sample,
        "gene": r["gene"],
        "copies": r["copies"],
        "estimate": r["estimate"],
        "confidence": r["confidence"],
        "status": r["status"],
        "ambiguous": r["ambiguous"],
        "method": r["method"],
        "lambda1": round(model.lambda1, 2),
        "mean_depth": call.support.get("mean_depth", ""),
        "alt_verdict": model.alt_verdict,
        "assembly": target_assembly,
        "n_pairs": stats.n_pairs,
        "notes": ";".join(filter(None, [r["notes"]] + stats.warnings)),
        "support": r["support"],
    }


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--inputs", type=Path, required=True,
                   help="list file of CRAM/BAM paths or URLs, optional index "
                        "and sample-name columns")
    p.add_argument("--reference", type=Path, required=True,
                   help="the FASTA the input CRAMs were compressed against, "
                        "needed to decode them; unused for BAM input")
    p.add_argument("--target", type=Path,
                   help="the reference to realign to, and the bwa index base. "
                        "Defaults to --reference. Point it at chm13v2.0.fa to "
                        "measure in T2T coordinates, where LILRA3 is on the "
                        "primary assembly")
    p.add_argument("--bwa-index", type=Path,
                   help="bwa index base, if it is not beside --target")
    p.add_argument("-o", "--output", type=Path, required=True,
                   help="LILRA6 copy number, TSV, one row per sample")
    p.add_argument("--genes",
                   help="comma-separated subset of the genes this target can "
                        "report, e.g. LILRA6. Defaults to all of them")
    p.add_argument("--threads", type=int, default=8,
                   help="threads per sample, for bwa and samtools (default 8)")
    p.add_argument("--jobs", type=int, default=1,
                   help="samples in parallel (default 1). Total load is "
                        "--jobs x --threads")
    p.add_argument("--outdir", type=Path, default=Path("results/lilra6"),
                   help="where per-sample QC lands")
    p.add_argument("--keep-realigned", action="store_true",
                   help="keep the realigned BAMs under <outdir>/realigned")
    args = p.parse_args()

    reference = args.reference.resolve()
    target = (args.target or args.reference).resolve()
    bwa_index = (args.bwa_index or target).resolve()

    if not reference.exists():
        raise SystemExit(f"{reference}: not found; run scripts/fetch_reference.sh")
    if not target.exists():
        raise SystemExit(f"{target}: not found; "
                         "run scripts/fetch_t2t_reference.sh for CHM13")

    # Which coordinate table the realigned BAM will be expressed in. Read from
    # the target's own .fai rather than inferred from its filename, and the
    # single most important thing to get right here: GRCh38 and CHM13 both call
    # the chromosome chr19, so measuring a CHM13 BAM with GRCh38 intervals does
    # not fail, it reads sequence ~3 Mb away and reports copy number for it.
    target_assembly, target_loci = realign.assembly_of_reference(target)
    genes = genes_for(target_assembly, args.genes)
    for ext in (".bwt", ".pac", ".sa", ".ann", ".amb"):
        if not Path(str(bwa_index) + ext).exists():
            raise SystemExit(
                f"{bwa_index}{ext}: not found — the bwa index is incomplete.\n"
                "Run scripts/fetch_bwa_index.sh")
    # Only GRCh38 has ALT contigs, so only GRCh38 needs the .alt file. On CHM13
    # its absence is correct rather than a misconfiguration, and warning about it
    # would train people to ignore the warning that matters.
    if target_assembly == "GRCh38" and not realign.index_is_alt_aware(bwa_index):
        # Not fatal — every call comes back an honest `not_measured` rather than
        # a wrong number. But it is the single mistake that turns a cohort into
        # apparent deletion homozygotes, and bwa reports it nowhere, so it is
        # said once, loudly, before any work happens.
        print(f"WARNING: {bwa_index}.alt is missing. The realignment will not "
              "be ALT-aware and every LILRA6 call will be `not_measured`.\n"
              "         Run scripts/fetch_bwa_index.sh.", file=sys.stderr)

    rows = read_inputs(args.inputs)
    args.outdir.mkdir(parents=True, exist_ok=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    print(f"{len(rows)} samples, {args.jobs} x {args.threads} threads, "
          f"realigning against {target_assembly} ({bwa_index}), "
          f"reporting {', '.join(genes)}", file=sys.stderr)

    results: list[dict] = []
    started = time.time()
    if args.jobs > 1:
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            futures = {
                pool.submit(call_one, r, str(reference), str(bwa_index),
                            args.outdir, args.threads, args.keep_realigned,
                            target_assembly, str(target), genes): r
                for r in rows
            }
            for i, fut in enumerate(as_completed(futures), 1):
                res = fut.result()
                results.append(res)
                _progress(i, len(rows), res)
    else:
        for i, r in enumerate(rows, 1):
            res = call_one(r, str(reference), str(bwa_index), args.outdir,
                           args.threads, args.keep_realigned, target_assembly,
                           str(target), genes)
            results.append(res)
            _progress(i, len(rows), res)

    results.sort(key=lambda r: r["sample"])
    all_rows = [row for res in results for row in res["rows"]]
    lilra6 = [r for r in all_rows if r["gene"] == "LILRA6"]

    _write_tsv(args.output, all_rows)
    (args.outdir / "qc.json").write_text(json.dumps(
        [{k: v for k, v in r.items() if k != "rows"} for r in results],
        indent=2, sort_keys=True))

    n_failed = sum(1 for r in results if not r["ok"])
    measured = [r for r in lilra6 if r["status"] == "measured"]
    print(f"\n{len(measured)}/{len(rows)} LILRA6 calls measured, "
          f"{n_failed} samples failed, {time.time() - started:.0f}s",
          file=sys.stderr)
    if measured:
        dist: dict = {}
        for r in measured:
            dist[r["copies"]] = dist.get(r["copies"], 0) + 1
        print("  LILRA6 CN: " + ", ".join(f"{k}:{dist[k]}" for k in sorted(dist)),
              file=sys.stderr)
    print(f"  {args.output}\n  {args.outdir / 'qc.json'}", file=sys.stderr)
    # A failed sample is a row in the output, not an exit code -- but a run where
    # everything failed is a configuration problem and should not look like a
    # success to a scheduler.
    return 1 if n_failed == len(rows) else 0


def _progress(i: int, total: int, res: dict) -> None:
    if res["ok"]:
        detail = "  ".join(f"{r['gene']}={r['copies']} ({r['status']})"
                           for r in res["rows"])
    else:
        detail = f"FAILED: {res['error'].splitlines()[0][:70]}"
    print(f"  [{i}/{total}] {res['sample']}: {detail} "
          f"[{res['elapsed_s']}s]", file=sys.stderr)


def _write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUTPUT_FIELDS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    raise SystemExit(main())
