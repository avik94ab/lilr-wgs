# `lilr_cn.py` — LILR copy number, one file, any cluster

Copy number for **LILRA3, LILRA6 and LILRB3** from aligned short-read WGS.

One Python file. One dependency: `samtools` on `PATH`. No Snakemake, no conda
environment, no scheduler, no reference panels, and nothing about the machine it
was written on. Copy it onto whatever cluster you have and run it.

```bash
curl -O https://raw.githubusercontent.com/avik94ab/lilr-wgs/main/portable/lilr_cn.py
chmod +x lilr_cn.py
./lilr_cn.py --check                       # is samtools here, and can it read URLs?
```

This is the copy-number half of [lilr-wgs](../README.md), extracted so it can be
used on its own. It produces the same numbers: `tests/test_portable.py` holds
every constant and every decision in lockstep with `src/lilrwgs/`, and on the
100-sample 1000 Genomes cohort the package was run against, the two agree on
**300/300 calls**, support values included.

It does *not* do the other half — phased allele sequences, the pangenome panels,
paralogue arbitration. Use the full pipeline for those.

## Usage

```bash
# a few files, 16 threads
lilr_cn.py -r GRCh38_full_analysis_set_plus_decoy_hla.fa -t 16 \
    -o cn_calls.tsv HG00096.final.cram HG00097.final.cram

# a cohort, from a manifest
lilr_cn.py -r GRCh38.fa -t 32 -m samples.tsv -o cn_calls.tsv --qc coverage.tsv

# straight from 1000 Genomes over HTTPS — htslib fetches only the containers
# that overlap the locus, so a sample costs ~13 MB against a 30 GB CRAM
lilr_cn.py -r GRCh38.fa -t 16 -m urls.tsv -o cn_calls.tsv
```

`-t/--threads` is the whole CPU story: give it what your job was allocated. With
several samples the budget runs them side by side; with one sample it all goes
into that sample's slice. `--workers` overrides the split if you want to.

A manifest is a TSV or CSV with `sample`/`sample_id` and `cram`/`bam`/`path`
columns — lilr-wgs's own manifests work unchanged — or a plain two-column file,
or just one path per line, in which case sample names come from the filenames.

Input can be CRAM or BAM, local or `https://`. A local, indexed BAM is read
directly; anything else is sliced first. `--reference` is required for CRAM and
must be **the FASTA the CRAM was compressed against**: a different GRCh38 will
decode most bases correctly and corrupt the rest, which is considerably worse
than failing.

## On a cluster

Nothing here knows what a scheduler is, so use whichever you have. Two shapes
work well.

**One job, several samples.** Simplest, and resumable — a rerun skips whatever is
already in `--keep-slices`.

```bash
#SBATCH -c 16
lilr_cn.py -r $REF -t $SLURM_CPUS_PER_TASK -m samples.tsv \
    --keep-slices $SCRATCH/slices -o cn_calls.tsv --qc coverage.tsv
```

**One array task per sample**, then concatenate.

```bash
#SBATCH --array=1-2504 -c 4
sed -n "${SLURM_ARRAY_TASK_ID}p" samples.tsv > $TMPDIR/one.tsv
lilr_cn.py -r $REF -t $SLURM_CPUS_PER_TASK -m $TMPDIR/one.tsv \
    -o out/cn.$SLURM_ARRAY_TASK_ID.tsv --quiet
```

SGE is the same with `-t 1-2504` and `$SGE_TASK_ID`; PBS with
`-J 1-2504` and `$PBS_ARRAY_INDEX`; with no scheduler at all,
`xargs -P` or GNU `parallel` over the manifest does the job.

**If your compute nodes have no outbound network** — common, and it does not
announce itself: a remote CRAM dies inside the job with `Destination address
required` while the same command works on the login node. Do the one networked
pass where there is a route, then call offline:

```bash
# login node
lilr_cn.py -r $REF -t 8 -m urls.tsv --slice-only --keep-slices $SCRATCH/slices
# compute nodes, no network needed, no reference needed
lilr_cn.py -t 32 -m urls.tsv --keep-slices $SCRATCH/slices -o cn_calls.tsv
```

Both invocations take the same manifest; the second finds the cached slices by
sample name and never touches the URLs.

**Cost.** A remote slice is ~13 MB and about a minute of wall clock, mostly
latency, so it parallelises well — 16 at a time sustains ~17 samples/minute.
Calling from a local slice is ~1–2 minutes of one thread, dominated by two
`samtools depth` passes over ~115 kb of control loci. Expect roughly 3 CPU-hours
per thousand samples, plus whatever the downloads cost.

## Output

One row per sample per gene, TSV:

| column | |
|---|---|
| `copies` | the call: an integer, or empty when nothing was measured |
| `estimate` | the continuous value before rounding |
| `confidence` | distance from the rounding boundary — 1.0 on an integer, 0.0 on a half |
| `method` | `unique_window_q20`, `alt_depth_q0` or `junction` |
| `status` | `measured` / `not_measured` / `failed` |
| `ambiguous` | the estimate landed within 0.2 of a half-integer |
| `notes` | why, when something is off |
| `support` | the raw counts the call was made from, as JSON |

`--qc` writes the per-sample coverage model: λ₁, the MAPQ-20 retention inside and
outside the alt placement, and `alt_verdict`. Read it. It is where you find out
whether the calls mean anything.

**`status` is not a boolean and `copies` is not a number you can default to 0.**
A failed measurement and a true zero are different values throughout:

- LILRA3 is deleted at ~24% allele frequency, so CN 0 is a common *true* state.
  A region query that failed and got reported as 0 manufactures deletions that
  look exactly like the real ones.
- If the alignment was not ALT-aware, every MAPQ-20 window in this cluster reads
  near zero — for everyone. Divided by a live baseline, a whole cohort comes out
  as confident LILRA6 deletion homozygotes. When that is detected, LILRA6 and
  LILRB3 are reported `not_measured`, `alt_verdict` says `not_alt_aware`, and no
  copy number is emitted. Do not fill it in.

## What it needs from your data

- **GRCh38, ALT-aware.** The 1000 Genomes 30× set, gnomAD/NYGC-style CRAMs, and
  anything else aligned with `bwa mem` *given the `.alt` file* qualify. Without
  it, LILRA3 is still callable (it is measured at MAPQ 0) but LILRA6 and LILRB3
  are refused rather than guessed. The script tells you which case you are in.
- **UCSC-style contig names** (`chr19`, not `19`).
- **The LRC alt contigs in the reference**, for LILRA3 by depth. If they are
  absent, LILRA3 falls back to the deletion-junction assay automatically and
  says so in `notes`.
- **Duplicate-marked, or not — either is fine.** Flagged duplicates are excluded;
  nothing is recomputed, since duplicate detection needs the whole library to
  judge and a 550 kb slice cannot.
- **Depth**: built and validated at 30×. Lower works with more ambiguous calls;
  λ₁ below 8 raises a warning.

## How the numbers are arrived at

- **λ₁**, the depth one haploid copy yields, is the median MAPQ-20 depth over
  four control loci in this same file, halved. Per sample and absolute — no
  cohort, no batch, one sample is enough.
- **LILRA6 and LILRB3** are ~97% identical over their 5' exons, so depth is taken
  only over the paralogue-unique 3' windows (2,900 bp and 1,881 bp), at MAPQ 20,
  and divided by λ₁.
- **LILRA3** is not on the GRCh38 primary assembly at all — the reference carries
  the deletion — so it is read as MAPQ-0 depth over the four alt contigs that do
  carry it, normalised on a MAPQ-0 baseline, *including supplementary alignments*
  (with `bwa mem -Y` that is where essentially all of its evidence sits). The
  deletion junction is measured independently as a check: they share no failure
  mode, and disagreement is reported rather than averaged away.
- **LILRA6 is cross-checked** against the pooled LILRA6+LILRB3 depth. That route
  is coarse, so disagreement adds a note and never changes the call.

## Validation

The package these numbers come from was scored against copy number derived from
HPRC assemblies for the **101 donors** who have both an assembly and a 1000
Genomes CRAM, and independently against
[JoGo-LILR](https://doi.org/10.1016/j.humimm.2025.111272) (Nagasaki et al., *Hum
Immunol* 2025), a published, methodologically unrelated LILRA6/LILRB3 caller, on
**200 further 1000 Genomes samples — 200/200 agreement on LILRA6**. See
[`../validation/README.md`](../validation/README.md) for the assembly comparison
and [`../validation/jogo_crosscheck.md`](../validation/jogo_crosscheck.md) for the
cross-check, including what the latter does not establish.

## Licence and citation

Same as the parent repository. If you use this, cite lilr-wgs and — if you lean
on the cross-validation — the JoGo-LILR paper above.
