#!/usr/bin/env python3
"""One sample's entire per-sample chain, in one process.

Slice the CRAM once, fit the coverage model, call copy number, recruit reads to
genes. Every intermediate stays on the compute node's local ``$TMPDIR``; only the
small outputs reach shared storage.

The fusion is inherited from `lilr-genotyper`, which learned it the hard way: as
separate Snakemake rules, the per-stage handoffs through a shared filesystem
stalled the scheduler on metadata latency — 90-second waits for outputs that
existed — and left the cluster idle. The intermediates here are a ~30 MB BAM and
11 per-gene BAMs, none of which anything downstream wants.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from lilrwgs import assign, cn, coverage, extract, loci  # noqa: E402
from lilrwgs.genotype import BOWTIE2_ARGS  # noqa: E402
from lilrwgs.shell import require  # noqa: E402


def align_to_panels(sample: str, r1: Path, r2: Path, panel_index: Path,
                    genes: list[str], work: Path, threads: int) -> dict[str, str]:
    """Align the extracted reads against each gene's pangenome panel.

    Recruitment, not calling. A 465-sequence panel means a divergent allele still
    attracts its own reads, which a single reference would lose — and losing them
    looks like low coverage at exactly the haplotypes that are most interesting.
    """
    require("bowtie2", "samtools")
    out: dict[str, str] = {}
    for gene in genes:
        bam = work / f"{gene}.bam"
        bt2 = subprocess.Popen(
            ["bowtie2", "-x", str(panel_index / gene), "-1", str(r1), "-2", str(r2),
             "-p", str(threads), *BOWTIE2_ARGS, "-S", "-"],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        proc = subprocess.run(["samtools", "view", "-b", "-o", str(bam), "-"],
                              stdin=bt2.stdout, capture_output=True)
        bt2.wait()
        if proc.returncode != 0 or not bam.exists():
            raise RuntimeError(f"{sample}/{gene}: panel alignment failed")
        # A killed or truncated conversion leaves a BAM that Snakemake would
        # accept as done and that would poison arbitration. A legitimately empty
        # alignment still has a valid header and passes.
        if subprocess.run(["samtools", "quickcheck", str(bam)]).returncode != 0:
            raise RuntimeError(f"{sample}/{gene}: BAM failed integrity check")
        out[gene] = str(bam)
    return out


def anchor_depths(sample: str, reads_dir: Path, locus_index: Path,
                  locus_refs: Path, work: Path, threads: int) -> dict[str, float]:
    """Median realigned depth at the copy-stable anchor genes.

    Measured on the same path the genotyper uses — recruited FASTQ realigned to
    the single per-locus reference — so the number is in the units a threshold
    will later be applied in. That is the whole point: lambda_1 from the CRAM and
    depth from a per-gene BAM are different quantities, and comparing them
    directly makes every gene look under-covered by the pipeline's own losses.
    """
    import statistics as stats

    out: dict[str, float] = {}
    for gene in coverage.EFFICIENCY_ANCHORS:
        r1 = reads_dir / f"{gene}_R1.fq.gz"
        r2 = reads_dir / f"{gene}_R2.fq.gz"
        if not r1.exists() or not r2.exists():
            continue
        bam = work / f"anchor_{gene}.bam"
        bt2 = subprocess.Popen(
            ["bowtie2", "-x", str(locus_index / gene), "-1", str(r1), "-2", str(r2),
             "-p", str(threads), *BOWTIE2_ARGS, "-S", "-"],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        subprocess.run(["samtools", "sort", "-o", str(bam), "-"],
                       stdin=bt2.stdout, capture_output=True)
        bt2.wait()
        if not bam.exists() or bam.stat().st_size == 0:
            continue
        res = subprocess.run(["samtools", "depth", "-a", str(bam)],
                             capture_output=True, text=True)
        depths = [int(ln.split("\t")[2]) for ln in res.stdout.splitlines()
                  if len(ln.split("\t")) >= 3]
        if depths:
            # Median, not mean: the per-locus references carry flanking sequence
            # that recruitment does not cover, and those zero-depth tails would
            # drag a mean down and inflate the apparent loss.
            out[gene] = float(stats.median(depths))
    return out


def write_gene_fastqs(filtered: dict[str, str], out_dir: Path,
                      threads: int) -> None:
    """Filtered per-gene BAM -> paired FASTQ, one pair per gene."""
    out_dir.mkdir(parents=True, exist_ok=True)
    for gene, bam in filtered.items():
        collated = Path(bam).with_suffix(".collate.bam")
        subprocess.run(["samtools", "collate", "-u", "-o", str(collated), bam],
                       capture_output=True, check=True)
        subprocess.run(
            ["samtools", "fastq", "-n", "-@", str(threads),
             "-1", str(out_dir / f"{gene}_R1.fq.gz"),
             "-2", str(out_dir / f"{gene}_R2.fq.gz"),
             "-0", "/dev/null", "-s", "/dev/null", str(collated)],
            capture_output=True, check=True)
        collated.unlink(missing_ok=True)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sample", required=True)
    p.add_argument("--cram", required=True)
    p.add_argument("--reference", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--panel-index", required=True)
    p.add_argument("--locus-index", required=True)
    p.add_argument("--locus-refs", required=True)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--buffer", type=int, default=assign.DEFAULT_BUFFER)
    p.add_argument("--shared-groups", default="LILRA6,LILRB3 LILRB1,LILRB4")
    p.add_argument("--keep-slice", default="0")
    p.add_argument("--genes-file")
    args = p.parse_args()

    out = Path(args.outdir)
    genes = ([ln.strip() for ln in open(args.genes_file) if ln.strip()]
             if args.genes_file else loci.GENES)
    groups = [tuple(g.split(",")) for g in args.shared_groups.split() if g]

    work = Path(os.environ.get("TMPDIR", "/tmp")) / f"lilrwgs_{args.sample}"
    if work.exists():
        shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    started = time.time()

    try:
        if args.cram.startswith(("http://", "https://", "s3://", "gs://")):
            extract.check_remote_support()

        # 1. The one remote pass.
        slice_bam = work / f"{args.sample}.slice.bam"
        slice_info = extract.slice_cram(args.sample, args.cram, args.reference,
                                        slice_bam, threads=args.threads,
                                        tmpdir=str(work))

        # 2. Coverage model, before anything depends on a depth number.
        model = coverage.measure(args.sample, str(slice_bam),
                                 reference=args.reference)
        cov_path = out / "coverage" / f"{args.sample}.json"
        cov_path.parent.mkdir(parents=True, exist_ok=True)
        cov_path.write_text(json.dumps({
            **model.as_row(),
            "gc_correction": model.gc_correction,
            "controls": [vars(c) for c in model.controls],
        }, indent=2))

        # 3. Copy number, per sample and absolute.
        cn_calls = cn.call_sample(args.sample, str(slice_bam), model,
                                  reference=args.reference)
        cn_path = out / "cn" / "per_sample" / f"{args.sample}.tsv"
        cn_path.parent.mkdir(parents=True, exist_ok=True)
        with cn_path.open("w", newline="") as fh:
            rows = [c.as_row() for c in cn_calls]
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows)

        # 4. LILR reads out, panels in, arbitration, per-gene FASTQ.
        r1, r2 = work / "lilr_R1.fq.gz", work / "lilr_R2.fq.gz"
        fq_stats = extract.to_fastq(args.sample, str(slice_bam), r1, r2,
                                    threads=args.threads, tmpdir=str(work),
                                    include_alts=slice_info["has_alt_contigs"])

        panel_bams = align_to_panels(args.sample, r1, r2, Path(args.panel_index),
                                     genes, work, args.threads)
        scores = {g: assign.paired_scores(b) for g, b in panel_bams.items()}
        assigned, astats = assign.arbitrate(scores, buffer=args.buffer,
                                            shared_groups=groups)
        filtered = assign.write_filtered_bams(panel_bams, assigned, work / "filtered")
        write_gene_fastqs(filtered, out / "reads" / args.sample, args.threads)
        n_shared = assign.write_shared_table(
            assigned, out / "reads" / args.sample / "shared_pairs.tsv")

        # 5. Calibrate the recruitment path, then rewrite the coverage model.
        # This has to come after the per-gene FASTQs exist, because it measures
        # the very path they came through.
        depths = anchor_depths(args.sample, out / "reads" / args.sample,
                               Path(args.locus_index), Path(args.locus_refs),
                               work, args.threads)
        model.efficiency, model.n_efficiency_anchors = \
            coverage.recruitment_efficiency(depths, model.lambda1)
        if model.n_efficiency_anchors == 0:
            model.warnings.append(
                "no anchor gene could be realigned, so the recruitment "
                "efficiency is unmeasured and left at 1.0; depth thresholds "
                "will be too strict by whatever the path actually loses")
        cov_path.write_text(json.dumps({
            **model.as_row(),
            "anchor_depths": depths,
            "gc_correction": model.gc_correction,
            "controls": [vars(c) for c in model.controls],
        }, indent=2))

        if args.keep_slice == "1":
            keep = out / "slices"
            keep.mkdir(parents=True, exist_ok=True)
            shutil.copy(slice_bam, keep / slice_bam.name)
            shutil.copy(str(slice_bam) + ".bai", keep / (slice_bam.name + ".bai"))

        stats_path = out / "qc" / "process" / f"{args.sample}.json"
        stats_path.parent.mkdir(parents=True, exist_ok=True)
        stats_path.write_text(json.dumps({
            "sample": args.sample,
            "slice": slice_info,
            "extraction": fq_stats.as_row(),
            "assignment": astats.as_row(),
            "n_shared_pairs": n_shared,
            "coverage": model.as_row(),
            "cn": [c.as_row() for c in cn_calls],
            "elapsed_s": round(time.time() - started, 1),
        }, indent=2))

        for warning in model.warnings + fq_stats.warnings:
            print(f"warning: {args.sample}: {warning}")
        print(f"{args.sample}: lambda1={model.lambda1:.1f} "
              f"eff={model.efficiency:.3f} "
              f"verdict={model.alt_verdict} "
              f"cn={{{', '.join(f'{c.gene}={c.copies}' for c in cn_calls)}}} "
              f"pairs={astats.n_total} dropped={astats.as_row()['dropped_fraction']:.3f} "
              f"in {time.time() - started:.0f}s")
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
