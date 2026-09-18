# CLAUDE.md

Guidance for Claude Code working in this repository.

## What this is

A Snakemake workflow that calls LILR copy number and phased allele sequences from
**1000 Genomes 30× srWGS CRAMs**. It is a rewrite of
[`lilr-genotyper`](https://github.com/avik94ab/lilr-genotyper), which does the
same job from targeted-capture FASTQs, and the rewrite exists because that
pipeline's depth constants are calibrated to several-hundred-× data. Read
`PLAN.md` before changing anything structural.

## Commands

```bash
micromamba create -y -f environment.yml && micromamba activate lilr-wgs
bash scripts/check_env.sh                    # asserts samtools has libcurl
bash scripts/fetch_reference.sh              # ~3.2 GB, once
bash scripts/fetch_bwa_index.sh              # ~5.3 GB, once; includes the .alt
bash scripts/build_indices.sh resources/gdna resources/gdna_index \
     resources/bundle/references resources/locus_index 8

python3 scripts/make_manifest.py --collection 2504 -o config/manifest.tsv
snakemake -s workflow/Snakefile --cores 32
snakemake -s workflow/Snakefile --profile profiles/sge     # Wynton

python -m pytest tests/ -q                   # 197 tests, no cluster needed

# The realigning front end: any GRCh38 CRAM, however it was aligned.
printf '%s\n' /data/*.final.cram > inputs.txt
python3 scripts/lilra6_cn.py --inputs inputs.txt \
    --reference resources/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    --threads 8 --jobs 4 -o lilra6_cn.tsv
qsub -t 1-25 scripts/lilra6_array.sh inputs.txt results/lilra6
```

Every stage is also a standalone CLI, so a failing sample can be debugged without
Snakemake:

```bash
export PYTHONPATH=src
python scripts/process_sample.py --sample HG00096 --cram <url|path> \
    --reference resources/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    --outdir results --panel-index resources/gdna_index \
    --locus-index resources/locus_index --locus-refs resources/bundle/references
python -m lilrwgs.genotype HG00096 LILRB1 results/reads/HG00096 \
    results/coverage/HG00096.json results/cn/per_sample/HG00096.tsv ...
```

## Architecture

```
CRAM ─► [process_sample] ─► coverage model, CN, per-gene reads
              └─► [genotype] ×11 ─► [summary]      [cohort_scale] (a leaf)

CRAM ─► [realign] ─► GRCh38 BAM ─► same coverage model, same CN ─► LILRA6 only
        (scripts/lilra6_cn.py — a second front end, not a second pipeline)
```

**There are two front ends and one everything-else.** `process_sample` reads
depth out of the alignment the CRAM arrived in. `scripts/lilra6_cn.py` extracts
the LRC and the control loci, converts them back to FASTQ, and realigns with
`bwa mem -Y` against the analysis set and its alt index. They differ only in how
the BAM is produced: `coverage.measure` and `cn.call_sample` take a BAM in
GRCh38 coordinates and do not care which aligner made it, which is exactly why
the two are comparable and why `validation/compare_realign.py` can compare them.
Do not fork the copy-number logic to serve the second path — if realigned input
needs different thresholds, the thresholds were wrong.

**What the realign path reports depends on the target assembly.** Against GRCh38
it is LILRA6 alone, via `cn.call_lilra6`; against CHM13 it is LILRA6 and LILRA3.
LILRB3 is measured wherever LILRA6 is, because LILRA6's only independent check is
the pooled LILRA6+LILRB3 depth, and is not reported. `genes_for()` decides this
from the target's `UNIQUE_WINDOWS`, and `--genes` narrows it.

LILRA3 is refused on GRCh38 rather than approximated: its depth route cannot
survive a regional extraction (§12) and its junction route scores 97/100 where
the CRAM-as-is path scores 100/100. A number this path calls less well than the
path beside it is worse than no number. If you add a gene here, say what it
scores against `cn_calls.tsv` first.

**There are two references and they do different jobs.** `--target` is what reads
are realigned to and measured in; `--reference` only *decodes* CRAM input, which
stores differences from a reference and cannot be read without one. BAM input
needs no `--reference`. Nothing after the slice touches it — mixing them up gives
you a GC curve, or a coordinate table, from the wrong assembly, and both produce
copy numbers rather than errors.

`src/lilrwgs/` is an importable package; `scripts/` holds drivers; `workflow/`
holds only the DAG. Pure logic (`depth_model`, the fitting functions in
`coverage`, `arbitrate` in `assign`) is separated from I/O deliberately so it can
be tested without a BAM.

`portable/lilr_cn.py` is a **deliberate duplicate** of the copy-number path —
one file, stdlib plus samtools, for people who want CN without the repository.
It re-declares every constant rather than importing one, so the failure mode is
drift: a threshold corrected in `src/lilrwgs/` and not there yields a script that
still runs and is quietly wrong. `tests/test_portable.py` compares the two sides
constant by constant and decision by decision; when it fails, the portable copy
is what is stale. Change both or neither.

**The DAG has no barrier.** `genotype` depends only on its own sample. This is
the payoff of measuring copy number against an absolute baseline, and it is worth
protecting: the predecessor's `cn_cohort` rule sits mid-DAG and forces a cohort to
run as one batch. If you find yourself adding a cohort-wide fit that something
downstream consumes, you have reintroduced it.

## Where the tests stop

Everything up to and including assignment is covered — `loci`, `loci_chm13`,
`coverage`, `depth_model`, `cn`, `realign`, `assign` and the portable duplicate,
250 tests. **`callability.py`, `genotype.py` and `sequences.py` have none**, which
is 1,023 lines carrying the callable track, variant calling and allele naming.

The 101-donor validation scored **copy number**. There is no evidence yet for any
claim about variant calls or allele assignments, and `PLAN.md`'s Phase 6 `[x]`
means the code was written, not that it was checked. Do not cite a sequence result
until `docs/variant_calling.md` §2 has been worked through; do not add a fixed
depth constant downstream of `depth_model` (§3 explains why `DP >= 6` is looser
than what is already enforced, not stricter).

The sharpest untested edge is `callability.gather_evidence()`: it is the junction
between arithmetic that is verified and I/O that is not, and its failure modes all
produce plausible depth with a wrong callable track.

## The depth model is the point

`src/lilrwgs/depth_model.py` is the reason this project exists. Three properties,
all load-bearing:

- **Per-sample.** Thresholds derive from λ₁, measured from control loci in the
  same CRAM. Never reintroduce a fixed DP constant — at 30× a `DP>=20` filter
  sits in the middle of the depth distribution and masks about half of every
  gene, which is what the predecessor does.
- **Ploidy-scaled.** A locus at *k* copies is the sum of *k* haploid negative
  binomials; NB is additive, so μ and *r* both scale by *k*.
- **Two-sided.** A position at 3× expected depth is an unseparated paralogue
  pile-up, and a "heterozygote" there is a paralogous sequence variant. The
  ceiling has no counterpart in the predecessor.

`effective_copies()` keeps the ceiling from backfiring on shared blocks: where
LILRA6/LILRB3 reads are kept in both genes, depth legitimately counts both, so the
expectation is the pair's combined copy number. Judged against the gene's own
copy number the gate would reject ~90% of the LILRA6 CDS in every sample.

**λ₁ has two units and they are not interchangeable.** It is measured on the CRAM
slice; it is *applied* to a realigned per-gene BAM, after extraction, panel
recruitment, arbitration and realignment have each dropped reads.
`recruitment_efficiency()` converts between them (0.813 on HG00096).
`CoverageModel.lambda_at()` applies it, so callers should use that rather than
`.lambda1` directly.

## Things that will bite you

- **`samtools view` with several regions emits a read once per region it
  overlaps.** This doubled control-locus depth, put λ₁ at 36.8 for a 30× library
  and halved every copy number — output that looked like a cohort of
  hemizygotes. Guarded twice: `loci.slice_intervals()` merges, and the view calls
  pass `-M`. Keep both.
- **Contig names ≠ file names.** `LILRB1_named.fa` contains a contig called
  `LILRB1`. Use `genotype._contig_name()`.
- **`samtools --version` prints two `Features:` lines.** libcurl is on the htslib
  one. Grep the whole output.
- **Relative paths break in pipeline stages**, which run in a scratch cwd.
  Resolve to absolute before handing anything to a tool.
- **`bwa` never says the `.alt` file is missing.** It looks for `{index}.alt`,
  and without it aligns without ALT-awareness — no warning, no non-zero exit.
  Every MAPQ-20 window in the LRC then reads near zero for every sample alike,
  which is a cohort of LILRA6 deletion homozygotes. `realign.index_is_alt_aware`,
  `scripts/check_env.sh` and `scripts/lilra6_array.sh` each check for the file
  because nothing downstream can tell that case from real data.
- **chr19 is called `chr19` in CHM13 and GRCh37 too**, and CHM13's LRC sits
  ~3 Mb right of GRCh38's — far enough to be a different LILR gene, close enough
  that every coordinate still resolves. `realign.source_assembly` and
  `realign.assembly_of_reference` both decide by chromosome 19's length
  (GRCh38 58,617,616 / CHM13 61,707,364), never by a filename.
- **chr19 is called `chr19` in GRCh37 too.** A name check passes and the
  extraction returns a different half-megabase. `realign.source_assembly` tells
  the two apart by chromosome 19's length (58,617,616 vs 59,128,983) and refuses
  rather than lifting over.
- **The module-provided samtools on Wynton has no libcurl.** Use the conda env.
- **Wynton compute nodes have no outbound internet.** A remote CRAM in a job dies
  with `Destination address required`; the login node is fine, so this is
  invisible until you submit. `scripts/stage_slices.py` does the one remote pass
  somewhere that has a route and rewrites the manifest to local slices. Then use
  `profiles/sge_staged`: `h_rt` picks the queue, and the Snakefile's 2 h default
  (sized for the remote read) buys nothing but `long.q`.
- **A shell block runs under `set -euo pipefail`.** `$PYTHONPATH` unset is the
  normal state on a compute node and fatal under `nounset`; every job in a cohort
  exited in under a second. Expand with `${PYTHONPATH:+:$PYTHONPATH}`.
  `tests/test_workflow.py` asserts this.
- **Killing Snakemake does not kill its jobs.** Orphans kept resubmitting a
  cached, pre-fix Snakefile while a new run fought them for the same outputs.
  Drain with `qstat -u $USER | awk 'NR>2{print $1}' | xargs -r qdel` first.
- **`--config` swallows positional targets.** Put `--jobs N` between them, and
  make targets absolute or they will not match an absolute `outdir`.

## Refusal is a feature

A failed measurement and a true zero are different values everywhere in this
codebase, and merging them is the most dangerous available bug because the output
stays plausible:

- LILRA3 CN 0 is a common true state (~24% deletion allele frequency), so a
  failed region query reported as 0 manufactures deletions.
- In a non-ALT-aware alignment every MAPQ-20 window in the LRC reads near zero.
  Divided by a live baseline that is a confident 0 copies, and a whole cohort of
  those looks exactly like LILRA6 deletion homozygotes.
  `CoverageModel.usable_mapq20` gates this; the status becomes `not_measured`.

Statuses are `measured` / `not_measured` / `failed`, and callability reasons are
`low_depth` / `high_depth` / `low_mapq` / `paralog_ambiguous` / `no_model`. Do not
collapse them into a boolean or a bare `N`.

## Paralogue handling

`assign.arbitrate()` sends each pair to its best gene by paired AS score. A tie
confined to one `shared_group` is kept in **every** tied gene, not dropped from
all of them: LILRA6/LILRB3 share a ~4.6 kb block with no gene-diagnostic 31-mers,
and a competitive discard deletes it from both at once (26% of the LILRB3 CDS in
the capture cohort). PING does the same for KIR2DL5A/B. Do not "fix" this.

The shared verdict is written to `shared_pairs.tsv` because reads are recruited
against a pangenome panel and called against a single reference, so there is a
FASTQ in between and BAM tags do not survive it. The predecessor lost the
information there.

## Two reference sets

| stage | reference | why |
|---|---|---|
| recruitment | `resources/gdna/{gene}.fasta` — 465-sequence HPRC panels | a divergent allele still attracts its own reads |
| calling | `resources/bundle/references/{gene}_named.fa` | variants need one coordinate system |

## Validation

`validation/build_truth.py` derives copy number and allele sequences from the
same HPRC panels the pipeline aligns against; **101 donors also have a 1000
Genomes CRAM**. Score with `validation/compare_cn.py`.

A donor in that overlap is aligned against a panel containing its own
haplotypes. `--leave-one-donor-out` writes donor-excluded panels; the rerun on
all 101 (`PLAN.md` §10) reproduced `cn_calls.tsv` byte for byte, because copy
number is measured on the CRAM slice at step 3 of `process_sample` and panels are
not opened until step 4. **Copy number cannot be inflated by panel circularity —
do not add a path that would make it so.** Recruitment is inflated, by up to 11%
at LILRB3, so the caveat is live for allele *sequence*, which has not been scored
yet. Report both runs, and say which claim each one supports.

## Conventions

- No large binaries in git. Panels (8.7 MB of FASTA) ship; indices and the 3.2 GB
  reference are built or fetched.
- Keep files host-agnostic; absolute Wynton paths belong in a site config.
- pytest from the start — the predecessor shipped with "no test runner and no
  linter" and the bugs listed above are the kind that finds.
- When a measurement contradicts a documented expectation, record the measurement
  in the code comment next to the constant it justifies, not only in a commit
  message.
