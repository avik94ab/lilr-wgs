#!/usr/bin/env python3
"""LILRA3, LILRA6 and LILRB3 copy number from aligned short-read WGS.

One file, one dependency (``samtools``), no cluster assumptions. Copy it onto
any machine that has samtools on PATH and run it:

    lilr_cn.py -r GRCh38_full_analysis_set_plus_decoy_hla.fa -t 16 \\
        -o cn_calls.tsv NA12878.cram HG00096.cram

It is the copy-number half of `lilr-wgs <https://github.com/avik94ab/lilr-wgs>`_,
extracted so that it can be used without Snakemake, without a conda environment
and without the 8.7 MB of pangenome panels the allele caller needs. The numbers
are the same ones: ``tests/test_portable.py`` asserts that every constant and
every decision here matches ``src/lilrwgs/``, and the two agree call-for-call on
the cohorts the package was validated against.

**Why this is not a depth-over-a-threshold script.** Three properties carry the
result, and dropping any one of them produces plausible, wrong output:

* Copy number is measured against **this sample's own** control loci, not
  against a cohort. λ₁ — the depth one haploid copy yields — comes from four
  housekeeping loci a few hundred kb away in the same file, so a single sample
  can be called on its own and no batch effect can move it.
* LILRA6 and LILRB3 are ~97% identical over their 5' exons, so they are read
  only over the **paralogue-unique windows** at their 3' ends, which exist only
  at MAPQ 20.
* MAPQ 20 in this cluster is conditional on the alignment having been
  **ALT-aware**. GRCh38 places nine LRC alt haplotypes at one primary interval
  covering every window used here; without the ``.alt`` file every MAPQ-20 count
  inside it reads near zero, for everyone. Divided by a live baseline that is a
  confident *nought copies* — a whole cohort of apparent LILRA6 deletion
  homozygotes. Three control loci sit outside that placement to detect it, and
  when it fires the answer is ``not_measured``, never 0.

A failed measurement and a true zero are different values in this output. LILRA3
is deleted at ~24% allele frequency, so CN 0 is a common true state and a failed
query reported as 0 manufactures deletions.

Requires: samtools >= 1.10 (>= 1.11 for `depth -G`), Python >= 3.9. For CRAMs
read over HTTPS the samtools build also needs libcurl; ``--check`` says so.
"""

from __future__ import annotations

import argparse
import concurrent.futures as futures
import csv
import json
import math
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

__version__ = "1.0.0"

# ---------------------------------------------------------------------------
# GRCh38 coordinates. Half-open, 0-based.
# ---------------------------------------------------------------------------

CHROM = "chr19"

# The slice taken out of each CRAM: LILRB3 to LILRB4 with ~65 kb of flank, which
# is far more than a 450 bp insert needs to keep pairs intact at the edges. It
# stops short of KIR3DL3 at 54,724,441 deliberately — the KIR cluster is the most
# copy-number-variable region in the genome and pulling it in buys nothing here.
LRC_SLICE = (CHROM, 54_150_000, 54_700_000)

# LILRA3 is *not on the primary assembly*: GRCh38's chr19 carries the common
# ~6.7 kb deletion, so the gene is annotated only on LRC alt contigs. Depth over
# LILRA3 is depth over these four intervals, and because they are near-identical
# those reads exist only at MAPQ 0 — the one deliberate exception to the MAPQ
# floor used everywhere else.
LILRA3_ALT = [
    ("chr19_GL949746v1_alt", 271_714, 278_472),
    ("chr19_GL949747v2_alt", 271_595, 278_354),
    ("chr19_GL949753v2_alt", 271_970, 278_728),
    ("chr19_KI270938v1_alt", 271_943, 278_701),
]

ALT_CONTIGS = [
    "chr19_GL949746v1_alt", "chr19_GL949747v2_alt", "chr19_GL949748v2_alt",
    "chr19_GL949749v2_alt", "chr19_GL949750v2_alt", "chr19_GL949751v2_alt",
    "chr19_GL949752v1_alt", "chr19_GL949753v2_alt", "chr19_KI270938v1_alt",
]

# Primary is the deleted allele, so a deleted chromosome's reads cross this base
# cleanly and a LILRA3-bearing one's soft-clip there. An assay with no failure
# mode in common with the depth route, which is the point of having it.
#
# This base was measured, not inherited: the value carried over from the
# predecessor, 54,296,977, is 28 bp left of where reads actually clip, and the
# assay duly returned ~0 for every donor regardless of copy number. Over 101 HPRC
# donors the pile-up sits at 54,297,005 for right-clipped reads and 54,297,010
# for left-clipped — 5 bp of microhomology at the Alu the breakpoint sits in.
LILRA3_JUNCTION = (CHROM, 54_297_005)
LILRA3_JUNCTION_MICROHOMOLOGY = 5

# Where an activating receptor and an inhibitory one genuinely differ: the 3'
# exons encoding the transmembrane and cytoplasmic regions. A whole-gene window
# at MAPQ 20 would measure mostly nothing.
UNIQUE_WINDOWS = {
    "LILRA6": [(CHROM, 54_236_589, 54_238_675), (CHROM, 54_239_049, 54_239_863)],
    "LILRB3": [(CHROM, 54_216_459, 54_217_051), (CHROM, 54_218_251, 54_219_540)],
}

# Gene bodies, for the pooled-pair cross-check only (UCSC ncbiRefSeqCurated).
GENE_SPANS = {
    "LILRA6": (CHROM, 54_236_589, 54_242_790),
    "LILRB3": (CHROM, 54_216_277, 54_223_007),
}

# (name, chrom, start, end, inside_alt_placement)
#
# Housekeeping loci within ~200 kb of the targets, none a known CNV, each < 2%
# ambiguous by a 100-mer self-similarity scan. The last three sit *outside* the
# primary placement of the LRC alt contigs, so no alt haplotype carries a second
# copy of them however bwa was run — which is what lets them answer whether the
# alignment was ALT-aware. LAIR1 is closer and was the obvious fifth control, but
# it is 28.4% ambiguous against LAIR2 and is deliberately not used.
CONTROLS = [
    ("PRPF31", CHROM, 54_115_410, 54_131_719, True),
    ("TMC4",   CHROM, 54_160_095, 54_173_250, True),
    ("MBOAT7", CHROM, 54_173_412, 54_189_882, True),
    ("TTYH1",  CHROM, 54_415_219, 54_436_904, True),
    ("PRKCG",  CHROM, 53_882_196, 53_907_652, False),
    ("CACNG6", CHROM, 53_991_148, 54_012_666, False),
    ("PPP6R1", CHROM, 55_229_778, 55_259_017, False),
]

FLANK = 1_000              # around each fetched interval; a 450 bp insert fits
MAPQ_STRICT = 20
MAPQ_ANY = 0
MIN_BASEQ = 13             # samtools' own mpileup default

# Secondary (0x100), QC-fail (0x200) and duplicate (0x400) are dropped at the
# slice. Supplementary (0x800) is deliberately KEPT: NYGC ran `bwa mem -Y`, which
# emits the ALT-contig hit of an ALT-aware alignment as a supplementary record,
# and LILRA3 is not on the primary assembly at all — on HG00099, 481 of 547
# records over one alt interval are supplementary. Excluding them made a
# two-copy donor read as a deletion homozygote.
EXCLUDE_FLAGS = 0x700
SUPPLEMENTARY = 0x800

# Verdict thresholds, from the LRC-vs-outside comparison.
MIN_Q20_LRC = 0.30         # below this, MAPQ-20 windows inside the LRC are dead
MAX_DILUTION = 1.50        # above this, reads are scattering across alt haplotypes

# Copy-number ranges seen in pangenome data. Estimates are clamped to these, so
# an estimate far outside becomes a low-confidence call at the boundary rather
# than an impossible integer.
CN_RANGE = {
    "LILRA3": (0, 2),      # biallelic insertion/deletion
    "LILRA6": (0, 6),
    "LILRB3": (0, 4),
}
AMBIGUOUS_BAND = 0.20

# A clip has to be long enough to be sequence rather than a trimmed base or two,
# and the aligned part long enough to place it.
MIN_CLIP = MIN_ANCHOR = 10
JUNCTION_TOLERANCE = 5
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")

REMOTE_PREFIXES = ("http://", "https://", "s3://", "gs://", "ftp://")


class ToolError(RuntimeError):
    """A tool exited non-zero, or a measurement could not be made."""


# ---------------------------------------------------------------------------
# Intervals
# ---------------------------------------------------------------------------


def merge_intervals(intervals):
    """Collapse overlapping and touching intervals, per contig.

    This is not tidiness. ``samtools view`` given several regions emits a read
    once *per region it overlaps*, so an interval contained in another silently
    doubles depth where they coincide — and three of the control loci sit inside
    the LRC slice. That inflates λ₁ by two and halves every copy number derived
    from it, which looks like a cohort of hemizygotes rather than like a bug.
    The ``-M`` flag below makes samtools emit each read once regardless; both
    guards are kept because they fail independently.
    """
    by_chrom = {}
    for chrom, start, end in intervals:
        by_chrom.setdefault(chrom, []).append((max(0, start), end))

    merged = []
    for chrom in sorted(by_chrom):
        spans = sorted(by_chrom[chrom])
        cur_start, cur_end = spans[0]
        for start, end in spans[1:]:
            if start <= cur_end:
                cur_end = max(cur_end, end)
            else:
                merged.append((chrom, cur_start, cur_end))
                cur_start, cur_end = start, end
        merged.append((chrom, cur_start, cur_end))
    return merged


def as_region(chrom, start, end):
    """0-based half-open -> a 1-based inclusive samtools region string."""
    return f"{chrom}:{start + 1}-{end}"


def slice_intervals(include_alts=True):
    """Everything any stage needs, merged, so the remote pass happens once.

    The alt-contig flank is generous because those interval boundaries come from
    41-mer anchoring rather than an annotation: being wrong by a few kb should
    cost coverage, not data.
    """
    intervals = [LRC_SLICE]
    if include_alts:
        intervals += [(c, max(0, s - 20_000), e + 20_000) for c, s, e in LILRA3_ALT]
    intervals += [(c, s - FLANK, e + FLANK) for _, c, s, e, _ in CONTROLS]
    chrom, pos = LILRA3_JUNCTION
    intervals.append((chrom, pos - 2_000, pos + 2_000))
    return merge_intervals(intervals)


# ---------------------------------------------------------------------------
# Running things
# ---------------------------------------------------------------------------


def run(cmd, *, text_input=None, check=True, env=None):
    """Run one command to completion with both streams captured.

    Failures carry the tool's own stderr: every failure mode here — a CRAM whose
    reference is wrong, a remote read that times out, a samtools too old for a
    flag — announces itself clearly in stderr and not at all in the exit code.
    """
    proc = subprocess.run([str(c) for c in cmd], capture_output=True, text=True,
                          input=text_input, env=env)
    if check and proc.returncode != 0:
        raise ToolError(f"{cmd[0]} exited {proc.returncode}\n"
                        f"  command: {' '.join(str(c) for c in cmd)}\n"
                        f"  stderr:  {proc.stderr.strip()[:1500]}")
    return proc


def cram_env(reference=None):
    """Environment for an htslib process reading a CRAM.

    ``REF_PATH`` is pinned rather than left to htslib's default, which points at
    the EBI CRAM reference registry. On a host that cannot reach it — or one that
    can, but slowly, a few thousand times — a silent outbound fetch per container
    hangs instead of failing, and the symptom is a job that merely looks slow.
    """
    env = dict(os.environ)
    env["REF_PATH"] = str(reference) if reference else "/dev/null/no_ref_registry"
    return env


def check_tools(samtools="samtools", *, remote=False):
    """Fail early and legibly, before a cohort discovers this the hard way."""
    if shutil.which(samtools) is None:
        raise ToolError(f"{samtools} is not on PATH. "
                        "conda install -c bioconda samtools, or module load samtools.")
    version = run([samtools, "--version"]).stdout
    if not remote:
        return version.splitlines()[0] if version else ""

    # `samtools --version` prints a Features line for samtools itself and another
    # for htslib; libcurl belongs to htslib. Checking only the first reports a
    # perfectly good installation as broken.
    features = [ln.strip() for ln in version.splitlines()
                if ln.strip().startswith("Features:")]
    if not any("libcurl=yes" in line for line in features):
        raise ToolError(
            f"{shutil.which(samtools)} is built without libcurl, so it cannot "
            "read CRAMs over HTTPS.\n"
            + "\n".join(f"  {line}" for line in features) +
            "\nUse a conda/bioconda build, or download the CRAMs first and pass "
            "local paths. Many HPC module builds omit libcurl."
        )
    return version.splitlines()[0] if version else ""


# ---------------------------------------------------------------------------
# Step 1: the one pass over the CRAM
# ---------------------------------------------------------------------------


def slice_cram(sample, source, reference, out_bam, *, threads=4,
               samtools="samtools", tmpdir=None):
    """Fetch every region any later step needs into one local, indexed BAM.

    htslib reads a remote CRAM through its ``.crai``, fetching only containers
    that overlap the requested regions, so a sample costs tens of megabytes
    against a CRAM of tens of gigabytes. The local BAM is also a resumption
    point: everything after it is cheap to redo, which is what ``--keep-slices``
    is for on a cluster whose compute nodes have no outbound network.
    """
    out_bam = Path(out_bam).resolve()
    out_bam.parent.mkdir(parents=True, exist_ok=True)
    # Stages may run in a scratch cwd, so every path handed to a tool is
    # absolute. A relative reference resolves against the scratch dir and fails
    # with "Failed to open reference file", which points at the reference rather
    # than at the cwd and is a confusing place to start looking.
    reference = Path(reference).resolve() if reference else None
    env = cram_env(reference)

    # Which contigs the CRAM actually has. Asking up front separates "aligned to
    # a primary-only reference" — a real, reportable fact — from "samtools
    # failed", which a bare error on an unknown region would conflate. Both
    # produce no LILRA3; only one of them is a bug.
    header_cmd = [samtools, "view", "-H"]
    if reference:
        header_cmd += ["-T", str(reference)]
    header = run(header_cmd + [source], env=env).stdout
    contigs = {part[3:] for ln in header.splitlines() if ln.startswith("@SQ")
               for part in ln.split("\t") if part.startswith("SN:")}

    has_alts = any(c in contigs for c in ALT_CONTIGS)
    regions = [as_region(*iv) for iv in slice_intervals(include_alts=has_alts)]
    regions = [r for r in regions if r.split(":")[0] in contigs]
    if not regions:
        raise ToolError(
            f"{sample}: none of the requested contigs are in the header.\n"
            f"  looked for: {CHROM}, {ALT_CONTIGS[0]}, ...\n"
            f"  header has: {', '.join(sorted(contigs)[:6])} ...\n"
            "This is usually Ensembl-style naming ('19' not 'chr19'); this "
            "script expects the UCSC-style names GRCh38 analysis sets use."
        )

    tmp = Path(tmpdir or os.environ.get("TMPDIR", tempfile.gettempdir()))
    tmp.mkdir(parents=True, exist_ok=True)
    half = max(1, threads // 2)
    with tempfile.TemporaryDirectory(prefix=f"lilrcn_{sample}_", dir=str(tmp)) as work:
        view = [samtools, "view", "-u", "-M", "-F", str(EXCLUDE_FLAGS),
                "-@", str(half)]
        if reference:
            view += ["-T", str(reference)]
        view += [source, *regions]
        sort = [samtools, "sort", "-@", str(half),
                "-T", str(Path(work) / "sort"), "-o", str(out_bam), "-"]

        # Wired explicitly rather than through a shell, because these arguments
        # include URLs and region strings and every one of them would become a
        # quoting problem. Both return codes are checked: a `view` killed by the
        # scheduler must not be masked by a `sort` that cheerfully writes a
        # valid, truncated, empty BAM.
        with open(Path(work) / "view.err", "w+") as verr, \
                open(Path(work) / "sort.err", "w+") as serr:
            p1 = subprocess.Popen(view, stdout=subprocess.PIPE, stderr=verr, env=env)
            p2 = subprocess.Popen(sort, stdin=p1.stdout, stderr=serr, env=env)
            p1.stdout.close()
            p2.wait()
            p1.wait()
            if p1.returncode or p2.returncode:
                verr.seek(0), serr.seek(0)
                raise ToolError(
                    f"{sample}: slicing failed (view={p1.returncode}, "
                    f"sort={p2.returncode})\n  {verr.read().strip()[:800]}\n"
                    f"  {serr.read().strip()[:800]}")

    run([samtools, "index", "-@", str(threads), str(out_bam)])
    n_records = int(run([samtools, "view", "-c", str(out_bam)]).stdout.strip() or 0)
    if n_records == 0:
        raise ToolError(
            f"{sample}: the slice is empty.\n  source: {source}\n"
            "The contig names matched, so this is most likely a reference "
            "mismatch on -T, or an index that does not cover these regions.")
    return {"bam": str(out_bam), "n_records": n_records,
            "n_regions": len(regions), "has_alt_contigs": has_alts}


# ---------------------------------------------------------------------------
# Step 2: what one haploid copy looks like in this sample
# ---------------------------------------------------------------------------


def _control_depths(bam, mapq, *, reference=None, samtools="samtools"):
    """Mean per-base depth over each control locus, keyed by name.

    ``-a`` matters: without it samtools omits zero-depth positions and a region
    that is half uncovered reports the mean of its covered half. That is the
    difference between "this control looks normal" and "half of this control is
    missing", which the model must not be blind to. Positions samtools never
    emits are counted as the zeros they are, so the denominator is always the
    locus span.

    ``-G 0x800`` drops supplementary records here, unlike at LILRA3: at unique
    control sequence a supplementary record is a misplaced fragment, not extra
    coverage, and counting it would inflate the baseline every threshold is
    expressed in.
    """
    bed = "".join(f"{c}\t{s}\t{e}\n" for _, c, s, e, _ in CONTROLS)
    cmd = [samtools, "depth", "-a", "-q", str(MIN_BASEQ), "-Q", str(mapq),
           "-G", "0x800", "-b", "/dev/stdin"]
    if reference:
        cmd += ["--reference", str(reference)]
    cmd.append(str(bam))
    out = run(cmd, text_input=bed).stdout

    totals = {name: 0 for name, *_ in CONTROLS}
    index = {}
    for name, chrom, start, end, _ in CONTROLS:
        index.setdefault(chrom, []).append((start, end, name))
    for line in out.splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        chrom, pos, depth = parts[0], int(parts[1]) - 1, int(parts[2])
        for start, end, name in index.get(chrom, ()):
            if start <= pos < end:
                totals[name] += depth
                break
    return {name: totals[name] / (end - start)
            for name, _c, start, end, _i in CONTROLS}


class CoverageModel:
    """λ₁, and whether to believe anything derived from it."""

    def __init__(self, sample):
        self.sample = sample
        self.lambda1 = 0.0
        self.lambda1_outside = 0.0
        self.q20_lrc = 0.0
        self.q20_outside = 0.0
        self.dilution = 0.0
        self.alt_verdict = "unknown"
        self.controls = []          # dicts: name, inside, mean_q0, mean_q20
        self.warnings = []

    @property
    def usable_mapq20(self):
        """Whether MAPQ-20 measurements inside the LRC mean anything.

        When they do not, the honest report is "not measured". Dividing a dead
        MAPQ-20 count by a live baseline produces a clean, confident zero, and a
        cohort of zeroes at LILRA6 looks exactly like a cohort of deletion
        homozygotes — a wrong answer that arrives looking like a right one.
        """
        return self.q20_lrc >= MIN_Q20_LRC

    def as_row(self):
        return {
            "sample": self.sample,
            "lambda1": round(self.lambda1, 3),
            "lambda1_outside": round(self.lambda1_outside, 3),
            "q20_lrc": round(self.q20_lrc, 4),
            "q20_outside": round(self.q20_outside, 4),
            "dilution": round(self.dilution, 4),
            "alt_verdict": self.alt_verdict,
            "usable_mapq20": self.usable_mapq20,
            "n_controls": len(self.controls),
            "warnings": ";".join(self.warnings),
        }


def alignment_verdict(q20_lrc, q20_outside, dilution):
    """Name what the diagnostics are saying, in words.

    - both retentions healthy, dilution ~1: an ALT-aware alignment, all good;
    - outside healthy, LRC collapsed: the alt contigs are eating MAPQ inside the
      placement — LILRA3 survives on the MAPQ-0 baseline, LILRA6 and LILRB3 do
      not and must be reported as not measured rather than as zero;
    - both collapsed: MAPQ is degraded genome-wide, a property of the sample or
      the aligner rather than of the alt contigs.
    """
    if q20_outside < MIN_Q20_LRC:
        return "mapq_degraded_everywhere"
    if q20_lrc < MIN_Q20_LRC:
        return "not_alt_aware"
    if dilution > MAX_DILUTION:
        return "alt_diluted"
    return "alt_aware"


def measure_coverage(sample, bam, *, reference=None, samtools="samtools"):
    """Fit the per-sample coverage model from the control loci.

    The controls sit at copy number 2, so λ₁ is half their depth. The median is
    taken across loci rather than the mean: one control overlapping an
    unannotated CNV in one sample should move the estimate by nothing, where with
    four loci a mean would move by a quarter of the error.
    """
    q0 = _control_depths(bam, MAPQ_ANY, reference=reference, samtools=samtools)
    q20 = _control_depths(bam, MAPQ_STRICT, reference=reference, samtools=samtools)

    model = CoverageModel(sample)
    for name, _chrom, _start, _end, inside in CONTROLS:
        model.controls.append({"name": name, "inside": inside,
                               "mean_q0": q0[name], "mean_q20": q20[name]})

    inside = [c for c in model.controls if c["inside"] and c["mean_q20"] > 0]
    outside = [c for c in model.controls if not c["inside"] and c["mean_q20"] > 0]
    if not inside and not outside:
        raise ToolError(
            f"{sample}: no control locus carried any depth in {bam}.\n"
            "If this is a pre-made slice, it was cut too narrowly — the controls "
            "sit deliberately outside the LILR cluster, from 53.88 to 55.26 Mb.")

    usable = inside or outside
    model.lambda1 = statistics.median(c["mean_q20"] for c in usable) / 2.0
    if outside:
        model.lambda1_outside = statistics.median(c["mean_q20"] for c in outside) / 2.0

    def retention(c):
        return c["mean_q20"] / c["mean_q0"] if c["mean_q0"] > 0 else 0.0

    model.q20_lrc = statistics.median([retention(c) for c in inside]) if inside else 0.0
    model.q20_outside = statistics.median([retention(c) for c in outside]) if outside else 0.0

    # Dilution: MAPQ-0 depth outside the alt placement over inside it. ~1 in an
    # ALT-aware alignment; well above 1 when reads that belong in the LRC are
    # being scattered across the alt haplotypes, which rescales every MAPQ-0
    # ratio normalised on an LRC control — LILRA3's, for instance.
    mean_in = statistics.median([c["mean_q0"] for c in inside]) if inside else 0.0
    mean_out = statistics.median([c["mean_q0"] for c in outside]) if outside else 0.0
    model.dilution = mean_out / mean_in if mean_in > 0 else 0.0
    model.alt_verdict = alignment_verdict(model.q20_lrc, model.q20_outside,
                                          model.dilution)

    if model.lambda1 <= 0:
        model.warnings.append("lambda1 is zero: no usable control depth")
    elif model.lambda1 < 8:
        model.warnings.append(
            f"lambda1 {model.lambda1:.1f} is low for a nominally 30x library; "
            "expect a high uncallable fraction")
    if model.alt_verdict == "not_alt_aware":
        model.warnings.append(
            "MAPQ 20 is dead inside the LRC but alive outside it: the alignment "
            "was not ALT-aware. LILRA6 and LILRB3 are not measurable by depth; "
            "they are reported as not measured rather than as zero")
    elif model.alt_verdict == "mapq_degraded_everywhere":
        model.warnings.append(
            "MAPQ 20 retention is low outside the alt placement too, so this is "
            "not an alt-contig problem — treat every MAPQ-filtered number here "
            "with suspicion")
    elif model.alt_verdict == "alt_diluted":
        model.warnings.append(
            f"dilution {model.dilution:.2f}: MAPQ-0 depth is markedly lower "
            "inside the alt placement than outside, so reads are scattering "
            "across the alt haplotypes")
    return model


# ---------------------------------------------------------------------------
# Step 3: copy number
# ---------------------------------------------------------------------------


class CNCall:
    """One gene's copy number in one sample, and how much to believe it."""

    def __init__(self, sample, gene, method=""):
        self.sample = sample
        self.gene = gene
        self.estimate = None       # continuous, before rounding
        self.copies = None
        self.confidence = 0.0      # 0 at a rounding boundary, 1 at an integer
        self.method = method
        self.status = "not_measured"   # measured | not_measured | failed
        self.support = {}
        self.notes = []

    @property
    def ambiguous(self):
        return self.status == "measured" and self.confidence < (1 - 2 * AMBIGUOUS_BAND)

    def as_row(self):
        return {
            "sample": self.sample,
            "gene": self.gene,
            "copies": "" if self.copies is None else self.copies,
            "estimate": "" if self.estimate is None else round(self.estimate, 3),
            "confidence": round(self.confidence, 3),
            "method": self.method,
            "status": self.status,
            "ambiguous": self.ambiguous,
            "notes": ";".join(self.notes),
            # The raw counts the call was made from. Carried because a copy
            # number without its evidence cannot be re-adjudicated later, and in
            # this cluster the calls that need re-adjudicating are exactly the
            # plausible-looking ones.
            "support": json.dumps(self.support, sort_keys=True),
        }


def integerise(estimate, gene):
    """Round an estimate to a copy number, with a confidence.

    Confidence is the distance from the nearest half-integer boundary, scaled so
    that landing on an integer gives 1.0 and landing on a boundary gives 0.0. It
    measures only how cleanly *this* estimate rounds — not whether λ₁ was right,
    whether the alignment was ALT-aware, or whether the gene was measurable at
    all. Those are separate fields precisely so a confident number derived from a
    broken measurement cannot masquerade as a good call.
    """
    lo, hi = CN_RANGE.get(gene, (0, 6))
    nearest = round(estimate)
    copies = int(min(hi, max(lo, nearest)))
    confidence = 1.0 - 2.0 * abs(estimate - nearest)
    # The clamping penalty applies only when clamping changed the answer: an
    # estimate of 2.047 at a gene capped at CN 2 rounds to 2 either way and is a
    # good call. Treating "outside the range" as "outside by any amount" scored
    # four of five correct LILRA3 calls at zero confidence, and since that flag
    # fires unevenly across copy-number classes it would bias any allele
    # frequency computed from the confident subset.
    if nearest != copies:
        confidence = 0.0
    return copies, max(0.0, confidence)


def _mean_depth(bam, intervals, mapq, *, reference=None, samtools="samtools",
                supplementary=False):
    """Mean per-base depth over a set of intervals, and the bases behind it.

    Returns None — never 0.0 — when the measurement could not be made. At LILRA3
    zero is the *expected* answer for a deletion homozygote at ~24% allele
    frequency, so conflating "no reads" with "could not ask" manufactures
    deletions out of failed queries.

    ``supplementary`` is False everywhere except LILRA3, where those records are
    the entire signal.
    """
    bed = "".join(f"{c}\t{s}\t{e}\n" for c, s, e in intervals)
    exclude = [] if supplementary else ["-G", "0x800"]
    cmd = [samtools, "depth", "-a", "-Q", str(mapq), *exclude, "-b", "/dev/stdin"]
    if reference:
        cmd += ["--reference", str(reference)]
    cmd.append(str(bam))
    try:
        out = run(cmd, text_input=bed).stdout
    except Exception:
        return None

    total = n = 0
    for line in out.splitlines():
        parts = line.split("\t")
        if len(parts) >= 3:
            total += int(parts[2])
            n += 1
    if n == 0:
        return None
    return total / n, n


def tally_junction(sam_lines, left, right):
    """Sort SAM records into breakpoint-clipped and breakpoint-spanning.

    The breakpoint has two ends, 5 bp apart across the microhomology: reads
    running into LILRA3 from the left flank end right-clipped at the first, reads
    coming back out of it start left-clipped at the second. Both are the bearing
    chromosome's signal, so both count as clipped. Soft-clipping elsewhere in a
    read is ordinary adapter or quality trimming and is ignored.
    """
    clipped = spanning = 0
    for line in sam_lines:
        fields = line.split("\t")
        if len(fields) < 6:
            continue
        start = int(fields[3]) - 1
        ops = CIGAR.findall(fields[5])
        if not ops:
            continue
        ref_len = sum(int(n) for n, op in ops if op in "MDN=X")
        end = start + ref_len
        if ref_len < MIN_ANCHOR:
            continue
        head, tail = ops[0], ops[-1]
        if (tail[1] == "S" and int(tail[0]) >= MIN_CLIP
                and abs(end - left) <= JUNCTION_TOLERANCE):
            clipped += 1
        elif (head[1] == "S" and int(head[0]) >= MIN_CLIP
                and abs(start - right) <= JUNCTION_TOLERANCE):
            clipped += 1
        elif start < left - MIN_ANCHOR and end > right + MIN_ANCHOR:
            spanning += 1
    return clipped, spanning


def junction_counts(bam, *, reference=None, samtools="samtools", window=300):
    """Reads clipped at the LILRA3 deletion junction, and reads spanning it."""
    chrom, left = LILRA3_JUNCTION
    right = left + LILRA3_JUNCTION_MICROHOMOLOGY
    region = f"{chrom}:{max(1, left - window)}-{right + window}"
    cmd = [samtools, "view", "-q", "1"]
    if reference:
        cmd += ["--reference", str(reference)]
    cmd += [str(bam), region]
    try:
        out = run(cmd).stdout
    except Exception:
        return None
    clipped, spanning = tally_junction(out.splitlines(), left, right)
    if clipped + spanning == 0:
        return None
    return clipped, spanning


def _call_unique_window(sample, gene, bam, model, *, reference, samtools):
    """LILRA6 or LILRB3: MAPQ-20 depth over the unique 3' windows, over λ₁."""
    call = CNCall(sample, gene, method="unique_window_q20")
    if model.lambda1 <= 0:
        call.status = "failed"
        call.notes.append("no coverage model")
        return call
    if not model.usable_mapq20:
        call.status = "not_measured"
        call.notes.append(
            f"MAPQ 20 is not usable in the LRC (q20_lrc={model.q20_lrc:.3f}, "
            f"verdict={model.alt_verdict}); reporting not measured rather than zero")
        return call

    measured = _mean_depth(bam, UNIQUE_WINDOWS[gene], MAPQ_STRICT,
                           reference=reference, samtools=samtools)
    if measured is None:
        call.status = "failed"
        call.notes.append("depth query failed")
        return call

    depth, n_bases = measured
    call.estimate = depth / model.lambda1
    call.copies, call.confidence = integerise(call.estimate, gene)
    call.status = "measured"
    call.support = {"mean_depth": round(depth, 2), "n_bases": n_bases,
                    "lambda1": round(model.lambda1, 2)}
    return call


def _call_lilra3(sample, bam, model, *, reference, samtools):
    """LILRA3 by two independent routes, preferring depth and checking it.

    Depth is preferred because it is a direct measurement over 6.7 kb, where the
    junction assay rests on clipping behaviour at a breakpoint inside an Alu. But
    the depth route needs the alt contigs in the file's reference, and where they
    are absent the junction is all there is — so the fallback is real, not
    decorative. The two share no failure mode, which is what makes their
    agreement worth more than either alone.
    """
    call = CNCall(sample, "LILRA3", method="alt_depth_q0")

    junction = junction_counts(bam, reference=reference, samtools=samtools)
    junction_estimate = None
    if junction is not None:
        clipped, spanning = junction
        # Over 101 HPRC donors this reads 0.00 at truth CN 0, a median 1.12 at
        # CN 1 and exactly 2.00 at CN 2. The 12% at CN 1 is the bearing
        # chromosome offering two breakpoints' worth of clipped reads against the
        # deleted one's single spanning window; left uncorrected, because it is
        # well inside the rounding band and a fitted fudge factor on a
        # cross-check would couple it to the route it exists to be independent of.
        junction_estimate = 2.0 * clipped / (clipped + spanning)
        call.support["junction_clipped"] = clipped
        call.support["junction_spanning"] = spanning
        call.support["junction_estimate"] = round(junction_estimate, 3)

    # MAPQ 0, normalised on a MAPQ-0 baseline: a ratio has to divide like by
    # like, and dividing a MAPQ-0 count by a MAPQ-20 baseline would inflate it by
    # whatever fraction of the baseline's reads are repeat-derived.
    lrc_q0 = [c["mean_q0"] for c in model.controls if c["inside"] and c["mean_q0"] > 0]
    measured = _mean_depth(bam, LILRA3_ALT, MAPQ_ANY, reference=reference,
                           samtools=samtools, supplementary=True)

    if measured is not None and lrc_q0:
        depth, n_bases = measured
        # Sum across the four intervals, divide by ONE interval's length. A read
        # from a LILRA3-bearing chromosome lands on exactly one of the four —
        # measured on HG00099, two of the contigs share zero read names out of
        # 518 and 464 — so the four counts partition the evidence rather than
        # replicating it, and the per-interval mean has to be multiplied back up.
        depth = depth * len(LILRA3_ALT)
        baseline_haploid = statistics.median(lrc_q0) / 2.0
        call.estimate = depth / baseline_haploid if baseline_haploid > 0 else None
        call.support["summed_depth"] = round(depth, 2)
        call.support["n_bases"] = n_bases
        call.support["n_alt_intervals"] = len(LILRA3_ALT)

    if call.estimate is None and junction_estimate is not None:
        call.estimate = junction_estimate
        call.method = "junction"
        call.notes.append(
            "no alt-contig depth (the file's reference has no LRC alt contigs); "
            "called from the deletion junction alone")

    if call.estimate is None:
        call.status = "failed"
        call.notes.append("neither the alt-contig depth nor the junction could be read")
        return call

    call.copies, call.confidence = integerise(call.estimate, "LILRA3")
    call.status = "measured"

    if junction_estimate is not None and call.method != "junction":
        disagreement = abs(junction_estimate - call.estimate)
        call.support["junction_disagreement"] = round(disagreement, 3)
        if disagreement > 0.75:
            call.notes.append(
                f"depth says {call.estimate:.2f} copies and the junction says "
                f"{junction_estimate:.2f}; these assays share no failure mode, so "
                "treat this call as unresolved")
            call.confidence = min(call.confidence, 0.3)
    return call


def _pair_check(calls, bam, model, *, reference, samtools):
    """Cross-check LILRA6 against the pooled LILRA6+LILRB3 depth.

    Coarse by construction: inverting the pooled ratio divides by the span weight
    and roughly doubles the noise. Read it for a systematic offset across a
    cohort rather than for per-sample agreement — which is why disagreement adds
    a note and never changes the call.
    """
    by_gene = {c.gene: c for c in calls}
    a6, b3 = by_gene.get("LILRA6"), by_gene.get("LILRB3")
    if not a6 or not b3 or a6.status != "measured" or b3.status != "measured":
        return
    if model.lambda1 <= 0:
        return

    gene_bodies = [GENE_SPANS[g] for g in ("LILRA6", "LILRB3")]
    spans = [e - s for _c, s, e in gene_bodies]
    measured = _mean_depth(bam, gene_bodies, MAPQ_ANY, reference=reference,
                           samtools=samtools)
    if measured is None:
        return
    depth, _ = measured

    # The span weight. Reads from both genes multi-map freely across both gene
    # bodies, so the pair's whole output spreads over the *sum* of the two spans
    # while one copy of one gene contributes coverage over one span. Dropping the
    # weight makes a 2+2 sample read as "pooled = 2", from which the implied
    # LILRA6 copy number comes out at zero for everybody.
    total_span = sum(spans)
    mean_span = total_span / len(spans)
    pooled = depth * total_span / (model.lambda1 * mean_span)
    implied_a6 = pooled - (b3.copies or 0)
    a6.support["pooled_pair_estimate"] = round(pooled, 2)
    a6.support["pair_implied_lilra6"] = round(implied_a6, 2)
    if a6.copies is not None and abs(implied_a6 - a6.copies) > 1.5:
        a6.notes.append(
            f"the pooled pair implies ~{implied_a6:.1f} LILRA6 copies against "
            f"{a6.copies} from the unique window; reads may be moving between "
            "the paralogues")


def call_sample(sample, bam, model, *, reference=None, samtools="samtools"):
    """Copy number for the three variable genes in one sample."""
    calls = [
        _call_unique_window(sample, "LILRA6", bam, model,
                            reference=reference, samtools=samtools),
        _call_unique_window(sample, "LILRB3", bam, model,
                            reference=reference, samtools=samtools),
        _call_lilra3(sample, bam, model, reference=reference, samtools=samtools),
    ]
    _pair_check(calls, bam, model, reference=reference, samtools=samtools)
    return calls


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

SAMPLE_SUFFIXES = (".final.cram", ".slice.bam", ".cram", ".bam", ".sam")


def sample_name(path):
    """A sample name from a path or URL, without the bioinformatics suffixes."""
    base = os.path.basename(str(path).split("?")[0].rstrip("/"))
    for suffix in SAMPLE_SUFFIXES:
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return os.path.splitext(base)[0]


def read_manifest(path):
    """(sample, source) pairs from a TSV or CSV.

    Accepts a header with ``sample``/``sample_id`` and ``cram``/``bam``/``path``
    columns — which is what lilr-wgs's own manifests look like — or a plain
    two-column file, or a bare list of paths.
    """
    rows = []
    delimiter = "," if str(path).endswith(".csv") else "\t"
    with open(path, newline="") as fh:
        lines = [ln for ln in fh.read().splitlines() if ln.strip()
                 and not ln.startswith("#")]
    if not lines:
        return rows

    header = lines[0].split(delimiter)
    lower = [h.strip().lower() for h in header]
    sample_col = next((i for i, h in enumerate(lower)
                       if h in ("sample", "sample_id", "id", "name")), None)
    path_col = next((i for i, h in enumerate(lower)
                     if h in ("cram", "bam", "path", "file", "url", "source")), None)
    if sample_col is not None and path_col is not None:
        body = lines[1:]
    else:
        sample_col, path_col, body = None, 0, lines

    for line in body:
        fields = line.split(delimiter)
        if len(fields) == 1:
            source = fields[0].strip()
            rows.append((sample_name(source), source))
            continue
        if sample_col is None:
            sample, source = fields[0].strip(), fields[1].strip()
        else:
            sample, source = fields[sample_col].strip(), fields[path_col].strip()
        if source:
            rows.append((sample or sample_name(source), source))
    return rows


def is_remote(source):
    return str(source).startswith(REMOTE_PREFIXES)


def has_index(source):
    p = Path(source)
    return any(p.with_suffix(p.suffix + ext).exists() or
               p.with_suffix(ext).exists()
               for ext in (".bai", ".csi", ".crai"))


def resolve_input(sample, source, args):
    """Where this sample's depth will actually be read from.

    ``local`` is an indexed BAM the region queries can read as it stands,
    ``cached`` a slice a previous run left in ``--keep-slices``, ``slice``
    anything that still has to be fetched.

    One function because the preflight checks and the worker have to agree. A
    cohort of URLs whose slices are all cached needs neither a reference nor a
    network, and refusing to start without them would send someone back to a
    login node for nothing — which is the whole staged workflow.
    """
    if not args.force_slice:
        if args.keep_slices:
            cached = Path(args.keep_slices) / f"{sample}.slice.bam"
            if cached.exists() and cached.stat().st_size > 0:
                return "cached", str(cached)
        if (not is_remote(source) and str(source).endswith((".bam", ".sam"))
                and has_index(source)):
            return "local", source
    return "slice", source


def process_one(sample, source, args, threads):
    """One sample, end to end: slice if needed, model, call. Never raises."""
    started = time.time()
    result = {"sample": sample, "source": source, "calls": [], "model": None,
              "error": None, "sliced": False}
    tmp_slice = None
    try:
        kind, slice_path = resolve_input(sample, source, args)
        if kind == "slice":
            if args.keep_slices:
                target = Path(args.keep_slices) / f"{sample}.slice.bam"
                target.parent.mkdir(parents=True, exist_ok=True)
            else:
                tmp_slice = tempfile.mkdtemp(prefix=f"lilrcn_{sample}_",
                                             dir=args.tmpdir or None)
                target = Path(tmp_slice) / f"{sample}.slice.bam"
            info = slice_cram(sample, source, args.reference, target,
                              threads=threads, samtools=args.samtools,
                              tmpdir=args.tmpdir)
            slice_path = info["bam"]
            result["sliced"] = True
            result["n_records"] = info["n_records"]

        if args.slice_only:
            result["elapsed_s"] = round(time.time() - started, 1)
            return result

        # The reference is only needed to decode CRAM; a BAM slice does not want
        # it, and passing one that does not match would be worse than passing
        # none at all.
        ref = args.reference if str(slice_path).endswith(".cram") else None
        model = measure_coverage(sample, slice_path, reference=ref,
                                 samtools=args.samtools)
        calls = call_sample(sample, slice_path, model, reference=ref,
                            samtools=args.samtools)
        result["model"] = model
        result["calls"] = calls
    except Exception as exc:                      # noqa: BLE001 - reported, not raised
        result["error"] = str(exc)
    finally:
        if tmp_slice:
            shutil.rmtree(tmp_slice, ignore_errors=True)
    result["elapsed_s"] = round(time.time() - started, 1)
    return result


def allocate(threads, n_samples):
    """Split one thread budget into (workers, samtools threads per worker).

    One knob, because one knob is what a user has: a CPU allocation. With several
    samples the budget goes to running them side by side, since each sample is a
    serial chain of short samtools calls that does not scale past a couple of
    threads; with one sample it all goes to that sample's slice.
    """
    threads = max(1, threads)
    workers = max(1, min(threads, n_samples))
    per_worker = max(1, threads // workers)
    return workers, per_worker


CN_FIELDS = ["sample", "gene", "copies", "estimate", "confidence", "method",
             "status", "ambiguous", "notes", "support"]
QC_FIELDS = ["sample", "lambda1", "lambda1_outside", "q20_lrc", "q20_outside",
             "dilution", "alt_verdict", "usable_mapq20", "n_controls", "warnings"]


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="lilr_cn.py",
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""examples:
  # a few CRAMs, 16 threads
  lilr_cn.py -r GRCh38.fa -t 16 -o cn.tsv A.cram B.cram C.cram

  # a cohort from a manifest (sample<TAB>path, header optional)
  lilr_cn.py -r GRCh38.fa -t 32 -m samples.tsv -o cn.tsv --qc qc.tsv

  # cluster with no outbound network on the compute nodes: fetch on the login
  # node, then call offline
  lilr_cn.py -r GRCh38.fa -t 8 -m urls.tsv --slice-only --keep-slices slices/
  lilr_cn.py -t 32 -m urls.tsv --keep-slices slices/ -o cn.tsv

  # one array task per sample (SLURM); $SLURM_ARRAY_TASK_ID is 1-based
  sed -n "${SLURM_ARRAY_TASK_ID}p" samples.tsv > $TMPDIR/one.tsv
  lilr_cn.py -r GRCh38.fa -t $SLURM_CPUS_PER_TASK -m $TMPDIR/one.tsv \\
      -o cn.$SLURM_ARRAY_TASK_ID.tsv
""")
    p.add_argument("inputs", nargs="*", metavar="CRAM|BAM|URL",
                   help="one or more aligned files; sample names come from the "
                        "filenames unless --manifest says otherwise")
    p.add_argument("-m", "--manifest",
                   help="TSV/CSV of samples: sample + cram/bam columns, a bare "
                        "two-column file, or one path per line")
    p.add_argument("-r", "--reference",
                   help="the FASTA the CRAM was compressed against. Required for "
                        "CRAM input; a *different* GRCh38 decodes most bases "
                        "correctly and corrupts the rest, which is worse than "
                        "failing. Not needed for BAM.")
    p.add_argument("-o", "--output", default="-",
                   help="copy-number TSV (default: stdout)")
    p.add_argument("--qc", help="per-sample coverage/QC TSV (optional but "
                                "recommended: it is where alt_verdict lives)")
    p.add_argument("-t", "--threads", type=int, default=4,
                   help="total CPU budget, split across samples (default: 4)")
    p.add_argument("--workers", type=int,
                   help="override how many samples run at once (default: "
                        "derived from --threads)")
    p.add_argument("--keep-slices", metavar="DIR",
                   help="write the ~13 MB per-sample slice here and reuse it on "
                        "a rerun; makes the run resumable")
    p.add_argument("--slice-only", action="store_true",
                   help="fetch slices and stop; pair with --keep-slices")
    p.add_argument("--force-slice", action="store_true",
                   help="re-slice even if a cached slice or local index exists")
    p.add_argument("--tmpdir", help="scratch for intermediates (default: $TMPDIR)")
    p.add_argument("--samtools", default="samtools", help="path to samtools")
    p.add_argument("--check", action="store_true",
                   help="verify samtools, including libcurl, and exit")
    p.add_argument("--quiet", action="store_true", help="no per-sample progress")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = p.parse_args(argv)

    if args.check:
        version = check_tools(args.samtools, remote=True)
        print(f"ok: {shutil.which(args.samtools)} ({version}), libcurl present")
        return 0

    samples = []
    if args.manifest:
        samples += read_manifest(args.manifest)
    samples += [(sample_name(s), s) for s in args.inputs]
    if not samples:
        p.error("no input: give one or more files, or --manifest")

    seen, unique = set(), []
    for sample, source in samples:
        if sample in seen:
            print(f"warning: duplicate sample {sample}; keeping the first",
                  file=sys.stderr)
            continue
        seen.add(sample)
        unique.append((sample, source))
    samples = unique

    if args.keep_slices:
        Path(args.keep_slices).mkdir(parents=True, exist_ok=True)

    # Judged on what will actually be opened, not on what the manifest says: the
    # second half of the staged workflow is handed a manifest of URLs whose
    # slices are already on disk, and it needs neither a reference nor a network.
    to_fetch = [source for sample, source in samples
                if resolve_input(sample, source, args)[0] == "slice"]
    if any(str(s).endswith(".cram") or is_remote(s) for s in to_fetch) \
            and not args.reference:
        p.error("CRAM input needs --reference (the FASTA it was compressed "
                "against). Slices already in --keep-slices do not.")
    if args.reference and not Path(args.reference).exists():
        p.error(f"reference not found: {args.reference}")
    if args.slice_only and not args.keep_slices:
        p.error("--slice-only without --keep-slices would fetch and discard")
    check_tools(args.samtools, remote=any(is_remote(s) for s in to_fetch))

    workers, per_worker = allocate(args.threads, len(samples))
    if args.workers:
        workers = max(1, args.workers)
        per_worker = max(1, args.threads // workers)
    if not args.quiet:
        print(f"lilr_cn {__version__}: {len(samples)} sample(s), {workers} at a "
              f"time, {per_worker} samtools thread(s) each", file=sys.stderr)

    started = time.time()
    results = []
    with futures.ThreadPoolExecutor(max_workers=workers) as pool:
        pending = {pool.submit(process_one, sample, source, args, per_worker):
                   sample for sample, source in samples}
        for done in futures.as_completed(pending):
            result = done.result()
            results.append(result)
            if args.quiet:
                continue
            n = len(results)
            if result["error"]:
                print(f"[{n}/{len(samples)}] {result['sample']}: FAILED "
                      f"{result['error'].splitlines()[0]}", file=sys.stderr)
            elif args.slice_only:
                print(f"[{n}/{len(samples)}] {result['sample']}: "
                      f"{'staged' if result['sliced'] else 'cached'} "
                      f"({result['elapsed_s']}s)", file=sys.stderr)
            else:
                model = result["model"]
                cn = ", ".join(f"{c.gene}={c.copies if c.copies is not None else '.'}"
                               for c in result["calls"])
                print(f"[{n}/{len(samples)}] {result['sample']}: "
                      f"lambda1={model.lambda1:.1f} {model.alt_verdict} "
                      f"cn={{{cn}}} ({result['elapsed_s']}s)", file=sys.stderr)
                for warning in model.warnings:
                    print(f"  warning: {result['sample']}: {warning}", file=sys.stderr)

    order = {sample: i for i, (sample, _) in enumerate(samples)}
    results.sort(key=lambda r: order[r["sample"]])

    if not args.slice_only:
        out = sys.stdout if args.output == "-" else open(args.output, "w", newline="")
        try:
            writer = csv.DictWriter(out, fieldnames=CN_FIELDS, delimiter="\t")
            writer.writeheader()
            for result in results:
                for call in result["calls"]:
                    writer.writerow(call.as_row())
        finally:
            if out is not sys.stdout:
                out.close()

        if args.qc:
            with open(args.qc, "w", newline="") as fh:
                writer = csv.DictWriter(fh, fieldnames=QC_FIELDS, delimiter="\t")
                writer.writeheader()
                for result in results:
                    if result["model"]:
                        writer.writerow(result["model"].as_row())

    failed = [r for r in results if r["error"]]
    if not args.quiet:
        n_called = sum(1 for r in results if r["calls"])
        print(f"done: {n_called} called, {len(failed)} failed, "
              f"{round(time.time() - started)}s", file=sys.stderr)
        for r in failed:
            print(f"  {r['sample']}: {r['error'].splitlines()[0]}", file=sys.stderr)
    # A cohort where every sample failed is a configuration error, not a result.
    return 1 if failed and len(failed) == len(results) else 0


if __name__ == "__main__":
    raise SystemExit(main())
