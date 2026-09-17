# lilr-wgs — project plan

Genotype the 11 LILR genes (copy number + phased allele sequences) from **1000 Genomes
Project 30× short-read WGS**, rather than from targeted capture.

This is a rewrite, not a port. The predecessor — `lilr-genotyper` — takes targeted-capture
FASTQs at several hundred ×, and almost every depth-related constant in it is calibrated to
that. At 30× WGS those constants are wrong in ways that do not announce themselves: a
`DP>=20` variant filter discards the majority of true heterozygous sites, and a
`depth < 15 → N` consensus mask blanks roughly half of every gene. The whole point of this
project is to replace fixed depth constants with a **per-sample depth model**, which srWGS
makes possible and capture data never did.

---

## 1. What changes, and why

| | `lilr-genotyper` (capture) | `lilr-wgs` (srWGS 30×) |
|---|---|---|
| input | `{sample}_R{1,2}.fastq.gz`, whole readset | GRCh38 CRAM, already aligned, ALT-aware |
| front end | 11 full-readset bowtie2 passes | one CRAM slice over the LRC → FASTQ |
| depth at LILR | hundreds ×, probe-driven, uneven | ~30×, near-uniform, GC-modulated |
| duplicates | must **not** be marked (probe geometry, not PCR) | already marked by NYGC; honour the flags |
| depth baseline | anchor LILR loci inside the same capture panel | genome-wide / chr19 control loci — a real, *external* baseline |
| CN calling | cohort GMM + unit fit; **needs ≥30 samples and a cohort barrier** | absolute, per-sample; cohort pass is a refinement, not a requirement |
| variant filter | `FMT/DP >= max(20, 10·CN)` — fixed | quantile of a fitted per-copy depth distribution, two-sided |
| consensus mask | `depth < max(20, 10·CN) → N` | callability track: floor **and** ceiling, MAPQ- and paralog-aware |
| validation | none in repo | 101 samples with HPRC assembly truth (see §6) |

Three of these deserve more than a table row.

### 1.1 Read depth is the subject, not a parameter

At 30× a diploid locus carries ~15 reads per haploid copy. Every threshold in the old
pipeline sits at or above that number, so it is not a filter any more — it is a mask. The
replacement is a model rather than a constant:

- **λ₁**, the expected depth per haploid copy, measured per sample from control loci that
  are outside the LILR cluster and outside GRCh38's LRC alt placement. GC-binned, because a
  30× library's depth varies by 20–30% across the GC range and the LILR genes are not at
  the genomic mean GC.
- At a locus of copy number *k*, expected depth is **λ_k = k·λ₁**, and observed depth is
  modelled as negative binomial around it (overdispersed Poisson; the dispersion is fitted
  per sample from the controls, not assumed).
- A position is **callable** when its depth falls inside a central interval of that
  distribution, and the interval is *two-sided*. This is the part the old pipeline had no
  concept of.

### 1.2 Excess depth is as diagnostic as missing depth

In a cluster of ~90%-identical paralogues, a position with 3× the expected depth is not
better evidence — it is a pile-up of reads from a paralogue that the assignment step failed
to separate, and a variant called there is a paralogous sequence variant masquerading as a
heterozygote. The old pipeline's one-sided floor accepted all of them silently. `lilr-wgs`
rejects depth above the upper quantile of λ_k and records *why* the position was dropped,
so a high-depth rejection and a low-depth rejection are never confused in the output.

### 1.3 The cohort barrier goes away

`lilr-genotyper` must run a cohort as one batch, because the per-haploid depth unit is fit
across all samples simultaneously (`cn_cohort` is a DAG barrier). With WGS, λ₁ is measurable
from each sample's own control loci, so copy number is callable from **one** sample. The
cohort pass survives as an optional refinement that sharpens the integer boundaries, and as
a QC instrument — but it never blocks the DAG. For 3,202 samples that matters.

---

## 2. What is reused, and from where

Not rewriting what already works:

| source | what | how it is used |
|---|---|---|
| `lilr-genotyper/workflow/scripts/filter_crossmapped.py` | paired-AS cross-map arbitration with `shared_groups` for the inseparable pairs (LILRA6/LILRB3, LILRB1/LILRB4) | ported nearly verbatim into `lilrwgs.assign`; the shared-group logic is load-bearing and stays |
| `lilr-genotyper/resources/gdna/*.fasta` | 11 HPRC pangenome gene panels, 231 donors | shipped unchanged; both the alignment target **and** the truth set (§6) |
| `lilr-genotyper` bundle `sequence_builder.py` | IUPAC-aware translation, CDS/cDNA extraction | ported into `lilrwgs.sequences` |
| `lilrCN_aou/resources/lilr_loci.json` + `lilr_cn.md` | GRCh38 coordinates: LILRA6/LILRB3 paralogue-unique windows, LILRA3 deletion junction at chr19:54,296,977, control loci inside and outside the alt placement, the ALT-awareness verdict | the coordinate basis of `lilrwgs.loci` and `lilrwgs.coverage` |
| **PING** (`Hollenbach-lab/PING`) | see below | |

From PING specifically:

- **`extractor_functions.R`** — one alignment pass against a combined reference, every later
  step reads the extracted FASTQs. `lilr-wgs` does the WGS equivalent: one CRAM slice.
- **`ping_copy.R` / `gc_functions.R`** — copy number as a *ratio to an internal reference
  locus* (KIR3DL3, fixed at CN 2), then thresholds to integers, with
  `manualCopyThresholds.csv` letting a human override a boundary the data leaves ambiguous.
  `lilr-wgs` adopts the override file wholesale; an auto-called boundary in a paralogous
  cluster should always be overridable, and PING is right to make that an artefact rather
  than a code edit.
- **The KFF probe idea** — presence/absence from exact k-mer counts in raw reads, normalised
  to the reference locus's probes. Alignment-free, nearly free to compute, and an
  independent check on LILRA3 presence that shares no failure mode with depth.
- **The precedent for not competitively filtering an inseparable pair** — PING comments out
  the negative filter for KIR2DL5A/B and shares one reference between them. That is the same
  decision `shared_groups` encodes for LILRA6/LILRB3.
- **`setup.minDP=8` / `final.minDP=20`** — a permissive threshold for building candidate
  genotypes and a strict one for finalising them. `lilr-wgs` keeps the two-stage shape but
  makes both stages relative to λ₁ instead of absolute.

---

## 3. Architecture

```
CRAM (GRCh38, ALT-aware)
   │
   ├─[1 extract]──► LRC read pairs ──► FASTQ ──┐
   │                                            │
   ├─[2 coverage]─► λ₁, dispersion, GC curve,   │
   │                ALT-aware verdict           │
   │                     │                      │
   └─[3 cn]──────────────┴──► CN per locus ─────┤   (per sample — no barrier)
        depth windows + junction + k-mer         │
                                                 ▼
                                      [4 assign] 11 panels + arbitration
                                                 │
                                                 ▼
                                      [5 genotype] HaplotypeCaller -ploidy CN
                                                 │   depth-model filters
                                                 │   whatshap phase/polyphase
                                                 ▼
                                      gDNA / CDS / cDNA / protein + callability track
                                                 │
                                                 ▼
                                      [6 report] summary + QC + cohort refinement
```

`src/lilrwgs/` is an importable package — every stage is a function with a CLI wrapper, so
it is unit-testable without Snakemake and debuggable without a cluster. The predecessor's
scripts were positional-argument CLIs with no tests; this is the main structural change.

---

## 4. Phases

Each phase ends with something that runs and something that is checked.

- [ ] **Phase 0 — skeleton.** Repo, plan, env spec, `.gitignore`, GitHub remote. Probe the
      Wynton toolchain (samtools with libcurl for remote CRAM, GATK, bowtie2, whatshap) and
      record what is actually available.
- [ ] **Phase 1 — resources and manifest.** `lilrwgs.loci` from the derived GRCh38
      coordinates. GRCh38 reference fetch script. Panel indices built, not shipped. Manifest
      builder turning the 1KGP sequence index into `sample_id → CRAM URL`.
- [ ] **Phase 2 — extraction.** CRAM → LRC read pairs → FASTQ, local path or EBI HTTPS.
      Verified end to end on HG00096/HG00097.
- [ ] **Phase 3 — the depth model.** ★ `lilrwgs.coverage` and `lilrwgs.depth_model`: λ₁,
      GC correction, dispersion, two-sided callability, the ALT-aware verdict. Unit tests
      against simulated depth where the answer is known.
- [ ] **Phase 4 — copy number.** Depth windows + LILRA3 junction + k-mer presence, combined
      into a per-sample integer call with an explicit confidence and an override file.
      Checked against HPRC truth (§6).
- [ ] **Phase 5 — per-gene assignment.** Panel alignment + ported cross-map arbitration,
      with the WGS-appropriate duplicate handling.
- [ ] **Phase 6 — genotyping.** HaplotypeCaller at ploidy = CN, depth-model filtering,
      phasing, consensus against the callability track, CDS/cDNA/protein.
- [ ] **Phase 7 — orchestration.** Snakemake DAG, SGE profile, resumability, the fused
      per-sample job pattern that the predecessor learned the hard way.
- [ ] **Phase 8 — validation.** The 101-sample HPRC comparison, CN and sequence, reported.
- [ ] **Phase 9 — documentation.** README, `CLAUDE.md`, method notes stating which claims
      rest on measurement and which on assumption.

---

## 5. Decisions taken up front

- **Snakemake**, matching the lab's existing workflow, with an SGE profile for Wynton.
- **A package plus thin CLIs**, not a pile of positional-argument scripts. pytest from
  Phase 0; the predecessor shipped with "no test runner and no linter" and it shows.
- **No large binaries in git.** The HPRC panels (8.7 MB of FASTA) ship; bowtie2 indices
  (29 MB) and the GRCh38 reference (3 GB) are built or fetched by `scripts/`.
- **Failure is never silently a zero.** A region query that fails and a region that is
  genuinely empty are different values with different status codes, all the way to the
  summary. At LILRA3, zero is a *common true answer* — the deletion is at ~15–20% allele
  frequency — so conflating the two would manufacture deletions.
- **Every masked position records its reason.** `low_depth`, `high_depth`, `low_mapq`,
  `paralog_ambiguous` are distinct in the callability track, never a bare `N`.

## 6. Validation design

The panel FASTA headers carry their donor: `>HG00097|hap1|LILRA6|copy1|HPRC|EUR|len=6358`.
231 HPRC donors are represented, and **101 of them also have a 1000 Genomes 30× CRAM**. For
those samples the HPRC assembly gives both truths at once:

- **copy number** — count the `copyN` entries per donor per gene across both haplotypes;
- **allele sequence** — the assembly contig itself, base for base.

So the pipeline can be run on the 1KGP CRAM for a sample whose answer is already known from
its assembly, and scored on both. The 101 samples are held out of nothing — they are not
training data, since no step is fitted on them — but any parameter chosen by looking at them
must be recorded as such in the method notes.

Care is needed in one place: the panels were built from those same assemblies, so a sample
in the overlap is aligning against a panel that contains its own haplotypes. That inflates
alignment rate and CN accuracy relative to an unseen sample. Phase 8 therefore scores twice —
once as-is, and once with the sample's own entries removed from the panel (leave-one-donor-out),
which is the number that generalises.

---

## 7. Open questions

Flagged rather than guessed; each is answered by a measurement in the phase noted.

1. **Is streaming CRAM slices from EBI fast enough for 2,504 samples, or must they be
   staged?** (Phase 2 — measure one sample, extrapolate.)
2. **Does the 30× depth support CN > 3 at LILRA6?** Capture data resolved CN 0–6. With
   λ₁ ≈ 15, separating CN 4 from CN 5 means separating 60× from 75× — feasible per position,
   marginal per sample. (Phase 4 — the HPRC truth set contains high-CN donors.)
3. **How much of LILRA6/LILRB3 is callable at 30× at all?** The paralogue-unique sequence is
   2,900 bp and 1,881 bp respectively; outside it the assignment step keeps shared-block
   reads in both genes, and whether that supports a haplotype call at this depth is
   unknown. (Phase 6/8.)
4. **Are the NYGC CRAMs ALT-aware in practice?** The functional-equivalence spec says yes;
   the tool measures it per sample rather than trusting it. (Phase 3.)
