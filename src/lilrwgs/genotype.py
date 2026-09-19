"""Variants, phasing and consensus for one (sample, gene).

Structurally this follows `lilr-genotyper`'s GATK pipeline, because that shape is
right: align to a single per-locus reference, call at ploidy = copy number, phase
read-backed, and emit one consensus per haplotype with the uncallable positions
masked. Three things differ, and all three are consequences of the input being
30x WGS rather than deep targeted capture.

**Duplicates are already gone.** The predecessor skips MarkDuplicates on purpose
— in targeted capture, pairs legitimately share start coordinates because probe
geometry drives where reads begin, and marking flagged 70-80% of them and crushed
depth to nothing. In WGS a shared start coordinate means what it usually means,
and the 1000 Genomes CRAMs arrive duplicate-marked by NYGC, so the flags are
honoured at the slice and nothing is recomputed here. Recomputing them on a
550 kb slice would be wrong anyway: duplicate detection needs the whole library
to judge.

**The variant filter is the depth model, not a constant.** The predecessor uses
``FMT/DP >= max(20, 10*CN)``. Here the callable positions come from
:mod:`lilrwgs.callability`, whose bounds vary per position with copy number,
local GC, MAPQ and whether the position sits in a shared paralogue block — and
are two-sided, so a pile-up is excluded as well as a dropout. The whole model
reaches bcftools as ``-T callable.bed`` without any of it having to survive being
written as a filter expression.

**Masking and filtering use one threshold, and it is the same one.** In the
predecessor these drifted apart at one point: the consensus mask sat at 15 and
the variant filter at 20, so a position with depth in [15, 20) survived the mask
but could never receive a variant call, and ``bcftools consensus`` emitted the
reference base — an unflagged reference-biased call rather than an honest N.
Deriving both from the same callability track makes that class of gap
unrepresentable.
"""

from __future__ import annotations

import json
import logging
import os
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

from . import callability
from .assign import DEFAULT_SHARED_GROUPS, read_shared_names
from .depth_model import HET_RATIO, MIN_DP_FINAL, MIN_DP_SETUP
from .sequences import extract_sequences
from .shell import ToolError, require, run

log = logging.getLogger(__name__)

GATK_JAVA_OPTS = ["--java-options", "-Xmx3g -Xms512m"]

# Wall-clock budget for phasing, per (sample, locus). `whatshap polyphase` cost
# climbs steeply with ploidy: measured on capture data at 2 threads, CN 3 took
# 69 min and CN 4 took 96 min. Unbounded, a handful of high-CN loci dominate a
# cohort's runtime, so phasing that overruns is abandoned and the locus falls
# back to the unphased consensus — recorded as phased_ok=False, not hidden.
PHASE_TIMEOUT = int(os.environ.get("LILRWGS_PHASE_TIMEOUT", 7200))

# bowtie2 settings, unchanged from the predecessor. Mixed and discordant pairs
# are deliberately NOT excluded: they recover paralogue-spanning reads at the
# LILRA1/LILRB2 D3 motif.
BOWTIE2_ARGS = [
    "-5", "3", "-3", "7", "-L", "20", "-i", "S,1,0.5",
    "--score-min", "L,0,-0.187", "-I", "75", "-X", "1000", "--no-unal",
]


def paralog_of(gene: str, groups=DEFAULT_SHARED_GROUPS) -> str | None:
    """The gene this one shares an inseparable block with, if any."""
    for members in groups:
        if gene in members and len(members) == 2:
            return members[0] if members[1] == gene else members[1]
    return None


@dataclass
class GenotypeResult:
    sample: str
    locus: str
    copies: int
    status: str = "ok"
    phased_ok: bool = False
    phasing_rate: float | None = None
    n_het: int = 0
    n_phased: int = 0
    callable_fraction: float = 0.0
    callability: dict = field(default_factory=dict)
    haplotypes: list[dict] = field(default_factory=list)
    error: str = ""

    def as_dict(self) -> dict:
        return {
            "sample": self.sample, "locus": self.locus, "cn": self.copies,
            "status": self.status, "phased_ok": self.phased_ok,
            "phasing_rate": self.phasing_rate,
            "n_het": self.n_het, "n_phased": self.n_phased,
            "callable_fraction": self.callable_fraction,
            **{f"cb_{k}": v for k, v in self.callability.items()
               if k.startswith("n_")},
            "haplotypes": self.haplotypes,
            "error": self.error,
        }


def _read_fasta(path: Path) -> str:
    return "".join(ln.strip() for ln in path.read_text().splitlines()
                   if not ln.startswith(">"))


def _contig_name(ref_fa: Path) -> str:
    """The contig name inside the reference, not the filename.

    These differ: the per-locus references are called ``LILRB1_named.fa`` and
    the sequence inside is ``LILRB1``. Deriving the name from the path puts the
    wrong contig in the callable BED, and GATK rejects the interval file with an
    error that names neither the BED nor the mismatch.
    """
    fai = ref_fa.with_suffix(ref_fa.suffix + ".fai")
    if fai.exists():
        first = fai.read_text().split("\n", 1)[0]
        if first.strip():
            return first.split("\t")[0]
    with ref_fa.open() as fh:
        for line in fh:
            if line.startswith(">"):
                return line[1:].split()[0]
    raise ToolError(f"{ref_fa}: no sequence header found")


def _ensure_ref_index(ref_fa: Path, gatk: str = "gatk") -> None:
    if not ref_fa.with_suffix(ref_fa.suffix + ".fai").exists():
        run(["samtools", "faidx", str(ref_fa)])
    if not ref_fa.with_suffix(".dict").exists():
        run([gatk, *GATK_JAVA_OPTS, "CreateSequenceDictionary", "-R", str(ref_fa)])


def process(
    sample: str,
    locus: str,
    reads_dir: str | os.PathLike,
    copies: int,
    paralog_copies: int,
    lambda1: float,
    dispersion: float,
    *,
    locus_ref_dir: str | os.PathLike,
    locus_index_dir: str | os.PathLike,
    protein_dir: str | os.PathLike,
    work_dir: str | os.PathLike,
    out_dir: str | os.PathLike,
    shared_table: str | os.PathLike | None = None,
    gc_lookup=None,
    alpha: float = 0.005,
    min_mapq_fraction: float = 0.5,
    het_ratio: float = HET_RATIO,
    min_dp: int | None = MIN_DP_FINAL,
    min_dp_setup: int | None = MIN_DP_SETUP,
    threads: int = 2,
    gatk: str = "gatk",
    whatshap: str = "whatshap",
    miniprot: str = "miniprot",
) -> GenotypeResult:
    """Genotype one gene in one sample."""
    require("bowtie2", "samtools", "bcftools")
    result = GenotypeResult(sample=sample, locus=locus, copies=copies)

    reads_dir = Path(reads_dir)
    r1, r2 = reads_dir / f"{locus}_R1.fq.gz", reads_dir / f"{locus}_R2.fq.gz"
    if not r1.exists() or not r2.exists():
        result.status = "skipped_no_reads"
        return result
    if copies == 0:
        # A true and common state at LILRA3, not a failure. Emitting nothing is
        # the honest output; emitting a reference-shaped consensus would invent
        # an allele the sample does not carry.
        result.status = "absent_cn0"
        return result

    work = Path(work_dir) / sample / locus
    work.mkdir(parents=True, exist_ok=True)
    out_sample = Path(out_dir) / sample
    out_sample.mkdir(parents=True, exist_ok=True)

    ref_fa = Path(locus_ref_dir) / f"{locus}_named.fa"
    _ensure_ref_index(ref_fa, gatk)
    ref_seq = _read_fasta(ref_fa)

    # 1. Align to the single per-locus reference.
    bam = work / "aligned.bam"
    bt2 = subprocess.Popen(
        ["bowtie2", "-x", str(Path(locus_index_dir) / locus),
         "-1", str(r1), "-2", str(r2), "-p", str(threads), *BOWTIE2_ARGS, "-S", "-"],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    subprocess.run(["samtools", "sort", "-@", str(threads), "-o", str(bam), "-"],
                   stdin=bt2.stdout, capture_output=True)
    bt2.wait()
    if not bam.exists() or bam.stat().st_size == 0:
        result.status = "no_alignment"
        return result

    rg_bam = work / "rg.bam"
    run([gatk, *GATK_JAVA_OPTS, "AddOrReplaceReadGroups",
         "-I", str(bam), "-O", str(rg_bam), "-RGID", sample, "-RGLB", sample,
         "-RGPL", "ILLUMINA", "-RGPU", sample, "-RGSM", sample,
         "--TMP_DIR", str(work)])
    run(["samtools", "index", str(rg_bam)])

    # 2. The depth model, per position.
    shared_names = read_shared_names(shared_table, locus) if shared_table else {}
    evidence = callability.gather_evidence(rg_bam, len(ref_seq),
                                           shared_names=shared_names)
    gc_by_pos = callability.gc_by_position(ref_seq) if gc_lookup else None
    # Two stages, after PING: permissive while discovering candidates, strict
    # at output. PING runs setup.minDP=8 / final.minDP=20; here it is 6 and 10,
    # lower because this is 30x WGS rather than capture at several hundred x.
    #
    # The predecessor also had two numbers -- masking the consensus at 15 and
    # filtering variants at 20 -- and that was a defect, because a position in
    # [15, 20) survived the mask and lost its call, so the consensus asserted a
    # base the VCF did not support. The difference here is which pair of things
    # share a threshold: discovery is looser than output, but the mask and the
    # variant filter are both driven by `calls_final`, so they cannot disagree.
    def _classify(min_dp):
        return callability.classify(
            evidence, copies=copies, lambda1=lambda1, dispersion=dispersion,
            paralog_copies=paralog_copies, alpha=alpha,
            min_mapq_fraction=min_mapq_fraction,
            gc_by_pos=gc_by_pos, gc_lookup=gc_lookup, min_dp=min_dp)

    calls_setup = _classify(min_dp_setup)
    calls = _classify(min_dp)          # final: the track, the mask, the filter

    contig = _contig_name(ref_fa)
    result.callability = callability.write_track(
        out_sample / f"{locus}.callability.tsv.gz".replace(".gz", ""),
        contig, evidence, calls)
    result.callable_fraction = result.callability["callable_fraction"]

    # What HaplotypeCaller may consider: the permissive set.
    bed = work / "callable.bed"
    n_callable = callability.callable_bed(bed, contig, calls_setup)
    # What survives into the VCF and the consensus: the strict set.
    final_bed = work / "final.bed"
    callability.callable_bed(final_bed, contig, calls)
    if n_callable == 0:
        result.status = "no_callable_positions"
        shutil.rmtree(work, ignore_errors=True)
        return result

    # 3. Call at ploidy = copy number, restricted to callable positions. Passing
    # -L here rather than filtering afterwards also saves HaplotypeCaller the
    # work of assembling regions that could never survive the filter.
    raw_vcf = work / "raw.vcf.gz"
    try:
        run([gatk, *GATK_JAVA_OPTS, "HaplotypeCaller",
             "-R", str(ref_fa), "-I", str(rg_bam), "-O", str(raw_vcf),
             "-ploidy", str(copies), "-L", str(bed), "--tmp-dir", str(work)])
    except ToolError as exc:
        # The callability track is already written and is the expensive part of
        # this function. Losing it because the caller failed would also lose the
        # evidence needed to work out why.
        result.status = "error"
        result.error = str(exc)
        return result

    # Indels are excluded from phasing and consensus, as upstream: they are not
    # reliable in a paralogous cluster at this depth and a wrong indel shifts
    # every downstream codon.
    filt_vcf = work / "filt.vcf.gz"
    # Allele balance, PING's `hetRatio`. A depth threshold cannot do this job:
    # at perfectly adequate depth, one or two reads from a 97%-identical
    # paralogue produce a heterozygote that passes any floor, and in this
    # cluster that is a more common error than thin coverage. The expression
    # requires the *minor* allele to carry at least `het_ratio` of the reads at
    # the site, so a 30x position with 2 alt reads (6.7%) is rejected while a
    # genuine het near 50% is kept.
    #
    # Applied only to sites that are actually called heterozygous -- a
    # homozygous-variant site has no minor allele to balance, and requiring one
    # would discard every one of them.
    ab = (f'(GT="het" & (FMT/AD[0:1])/(FMT/DP) >= {het_ratio}) '
          f'| GT!="het"')
    run(["bcftools", "view", "-V", "indels", "-T", str(final_bed),
         "-i", ab,
         "-Oz", "-o", str(filt_vcf), str(raw_vcf)])
    run(["bcftools", "index", "-f", str(filt_vcf)])

    # 4. Phase. `whatshap phase` is diploid-only; above CN 2 it silently produces
    # an unphased VCF that bcftools consensus then splits arbitrarily, so
    # polyphase is the correct entry point there.
    consensus_vcf = filt_vcf
    if copies >= 2:
        phased = work / "phased.vcf.gz"
        cmd = ([whatshap, "phase"] if copies == 2 else
               [whatshap, "polyphase", "--ploidy", str(copies),
                "--threads", str(threads)])
        try:
            run(cmd + ["-o", str(phased), "--reference", str(ref_fa),
                       "--ignore-read-groups", str(filt_vcf), str(rg_bam)],
                timeout=PHASE_TIMEOUT)
            run(["bcftools", "index", "-f", str(phased)])
            consensus_vcf = phased
            result.phased_ok = True
        except subprocess.TimeoutExpired:
            log.warning("%s/%s: phasing exceeded %ss at cn=%s; using unphased",
                        sample, locus, PHASE_TIMEOUT, copies)
        except Exception as exc:                       # noqa: BLE001
            log.warning("%s/%s: phasing failed (%s); using unphased",
                        sample, locus, exc)

    result.n_het, result.n_phased = _phasing_counts(consensus_vcf)
    result.phasing_rate = (result.n_phased / result.n_het
                           if result.n_het else None)

    # 5. One consensus per haplotype, masked by the same track that bounded the
    # variant call — so a position cannot be unfilterable and unmasked at once.
    for hap in range(1, copies + 1):
        hap_fa = work / f"hap{hap}.fa"
        run(["bcftools", "consensus", "-f", str(ref_fa), "-H", str(hap),
             "-o", str(hap_fa), str(consensus_vcf)])
        gdna = callability.mask_sequence(_read_fasta(hap_fa), calls)
        result.haplotypes.append(
            _finish_haplotype(sample, locus, hap, gdna, work, out_sample,
                              Path(protein_dir), miniprot))

    shutil.rmtree(work, ignore_errors=True)
    return result


def _phasing_counts(vcf: Path) -> tuple[int, int]:
    try:
        out = run(["bcftools", "view", str(vcf)]).stdout
    except ToolError:
        return 0, 0
    n_het = n_phased = 0
    for line in out.splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 10:
            continue
        gt = fields[9].split(":")[0]
        alleles = gt.replace("|", "/").split("/")
        if len(set(alleles)) > 1:
            n_het += 1
            if "|" in gt:
                n_phased += 1
    return n_het, n_phased


def _finish_haplotype(sample: str, locus: str, hap: int, gdna: str,
                      work: Path, out_sample: Path, protein_dir: Path,
                      miniprot: str) -> dict:
    """miniprot for CDS coordinates, then slice and translate."""
    hap_fa = work / f"hap{hap}_masked.fa"
    hap_fa.write_text(f">{sample}_{locus}_hap{hap}\n{gdna}\n")

    gdna_cds = cdna = protein = ""
    try:
        res = subprocess.run([miniprot, "--gff", str(hap_fa),
                              str(protein_dir / f"{locus}_protein.faa")],
                             capture_output=True, text=True, timeout=180)
        coords = _parse_gff(res.stdout)
        if coords:
            gdna_cds, cdna, protein = extract_sequences(gdna, coords)
    except Exception as exc:                            # noqa: BLE001
        log.warning("%s/%s/hap%s: miniprot failed (%s)", sample, locus, hap, exc)

    for kind, seq in (("gdna", gdna), ("gdna_cds", gdna_cds),
                      ("cdna", cdna), ("protein", protein)):
        if seq:
            (out_sample / f"{locus}_hap{hap}_{kind}.fasta").write_text(
                f">{sample}|{locus}|hap{hap}\n{seq}\n")

    return {
        "hap": hap,
        "gdna_len": len(gdna), "gdna_cds_len": len(gdna_cds),
        "cdna_len": len(cdna), "protein_len": len(protein),
        "gdna_n_pct": round(100 * gdna.count("N") / max(len(gdna), 1), 2),
        "protein_x_pct": round(100 * protein.count("X") / max(len(protein), 1), 2),
    }


def _parse_gff(text: str) -> dict | None:
    """mRNA span, strand and CDS exons from miniprot GFF, 0-based half-open."""
    mrna, strand, cds = None, "+", []
    for line in text.splitlines():
        if line.startswith("#"):
            continue
        parts = line.rstrip().split("\t")
        if len(parts) < 7:
            continue
        start, end = int(parts[3]) - 1, int(parts[4])
        if parts[2] == "mRNA":
            mrna, strand = (start, end), parts[6]
        elif parts[2] == "CDS":
            cds.append((start, end))
    if mrna is None:
        return None
    cds.sort()
    return {"mrna_start": mrna[0], "mrna_end": mrna[1],
            "strand": strand, "exons": cds}


def main(argv: list[str] | None = None) -> int:
    import argparse

    p = argparse.ArgumentParser(description="Genotype one (sample, gene)")
    p.add_argument("sample")
    p.add_argument("locus")
    p.add_argument("reads_dir")
    p.add_argument("coverage_json")
    p.add_argument("cn_tsv")
    p.add_argument("--locus-refs", required=True)
    p.add_argument("--locus-index", required=True)
    p.add_argument("--protein-dir", required=True)
    p.add_argument("--work-dir", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--shared-table")
    p.add_argument("--marker", required=True)
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--alpha", type=float, default=0.005)
    args = p.parse_args(argv)

    cov = json.loads(Path(args.coverage_json).read_text())
    dispersion = float("inf") if cov["dispersion"] == "inf" else float(cov["dispersion"])
    gc = {int(k): float(v) for k, v in (cov.get("gc_correction") or {}).items()}

    # lambda_1 is measured on the CRAM slice; depth here is measured on a
    # realigned per-gene BAM, after extraction, panel recruitment, arbitration
    # and realignment have each dropped reads. The efficiency factor converts
    # between the two, and without it every gene reads as under-covered by the
    # pipeline's own losses -- which is the predecessor's failure mode arriving
    # by a different route.
    efficiency = float(cov.get("efficiency", 1.0))
    lambda1 = float(cov["lambda1"]) * efficiency

    def gc_lookup(value: float) -> float:
        from .coverage import gc_bin
        return lambda1 * gc.get(gc_bin(value), 1.0)

    copies, paralog_copies = _copies_from(args.cn_tsv, args.sample, args.locus)

    try:
        result = process(
            args.sample, args.locus, args.reads_dir, copies, paralog_copies,
            lambda1, dispersion,
            locus_ref_dir=args.locus_refs, locus_index_dir=args.locus_index,
            protein_dir=args.protein_dir, work_dir=args.work_dir,
            out_dir=args.out_dir, shared_table=args.shared_table,
            gc_lookup=gc_lookup if gc else None, alpha=args.alpha,
            threads=args.threads,
        )
    except Exception as exc:                            # noqa: BLE001
        result = GenotypeResult(sample=args.sample, locus=args.locus,
                                copies=copies, status="error", error=str(exc))

    Path(args.marker).parent.mkdir(parents=True, exist_ok=True)
    Path(args.marker).write_text(json.dumps(result.as_dict(), indent=2))
    print(f"{args.sample}/{args.locus}: status={result.status} cn={copies} "
          f"callable={result.callable_fraction:.3f}")
    return 0


def _copies_from(cn_tsv: str, sample: str, locus: str) -> tuple[int, int]:
    """Copy number for this locus, and for its inseparable partner.

    Genes not called from data are fixed at 2. A gene whose call is
    ``not_measured`` also falls back to 2 rather than to 0 — reporting an
    unmeasurable gene as absent would turn a measurement failure into a
    biological claim.
    """
    import csv

    calls: dict[str, dict] = {}
    with open(cn_tsv) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["sample"] == sample:
                calls[row["gene"]] = row

    def copies_of(gene: str) -> int:
        row = calls.get(gene)
        if not row or row.get("status") != "measured" or row.get("copies") in ("", None):
            return 2
        return int(row["copies"])

    partner = paralog_of(locus)
    return copies_of(locus), copies_of(partner) if partner else 0


if __name__ == "__main__":
    raise SystemExit(main())
