#!/usr/bin/env python3
"""Derive the CHM13v2.0 intervals in `lilrwgs.loci_chm13`, and show the working.

    python3 scripts/derive_loci_chm13.py \
        --reference resources/reference/chm13v2.0.fa \
        --gff resources/reference/chm13v2.0_RefSeq_Liftoff_v5.1.gff3.gz \
        -o resources/chm13_loci.tsv

Nothing in `loci_chm13.py` is lifted over from GRCh38. Liftover is least
trustworthy exactly where this pipeline works — segmental duplications with
97%-identical paralogues — and a coordinate that lands one gene to the left
produces copy numbers rather than errors. Everything is derived from the CHM13
sequence and annotation directly, by this script, so it can be rerun and diffed.

Three kinds of interval, three different derivations:

**Gene bodies and control loci** come from the RefSeq Liftoff annotation, taking
the span per gene symbol. Note that LILRA3 is *not* in that file: it is lifted
from the GRCh38 primary assembly, which carries the deletion and so has no
LILRA3 to lift. Its absence there is an artefact of the annotation's provenance
and says nothing about the sequence.

**LILRA3** is located by aligning the calling reference to the region instead.

**The paralogue-unique windows** are *transferred* from GRCh38 by aligning their
sequence, and then checked. They are deliberately not re-derived from scratch,
and the reason is a measurement: running the documented criterion — a 100-mer is
ambiguous if another within 3 mismatches exists anywhere in the LRC on either
strand — against GRCh38 reproduces the *starts* of `loci.UNIQUE_WINDOWS` exactly
(54,236,589 and 54,218,251) but runs longer and finds extra windows, 4,003 bp
against 2,900 at LILRA6 and 3,539 against 1,881 at LILRB3. So `lilrCN_aou`
applied something stricter than the criterion as written, and LILRA6's validated
accuracy rests on the stricter version, not on this one.

Transferring the sequence keeps the definition that was validated and changes
only the assembly it is expressed in. The independent scan is still run, as a
check rather than as the source: a transferred window that is *not* unique in
CHM13 is reported loudly, because that would mean the two assemblies differ in a
way that matters to the measurement.

The scan is exact rather than heuristic. With at most 3 mismatches in 100 bp,
the pigeonhole principle puts at least one exact 25-mer in any twin, so indexing
every 25-mer of the region and verifying the candidates it proposes cannot miss
a twin that a full quadratic scan would find.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

# The LRC and its neighbourhood in CHM13, the search space a 100-mer is asked to
# be unique within. Chosen to match the 1.2 Mb GRCh38 window `lilrCN_aou` used
# (chr19:53.9-55.1 Mb), which spans PRKCG through PPP6R1 -- the same genes bound
# it here, so the question being asked of each 100-mer is the same question.
LRC_SEARCH = ("chr19", 56_980_000, 58_180_000)

WINDOW = 100        # the 100-mer of the criterion
MAX_MISMATCH = 3    # "a twin within 3 mismatches"
# floor(WINDOW / (MAX_MISMATCH + 1)): the longest stretch a twin is guaranteed to
# share exactly. Making this larger would miss twins; smaller only costs time.
SEED = WINDOW // (MAX_MISMATCH + 1)

# Genes whose spans are taken from the annotation.
LILR_GENES = ["LILRB3", "LILRA6", "LILRB5", "LILRB2", "LILRA5",
              "LILRA4", "LILRA2", "LILRA1", "LILRB1", "LILRB4"]
CONTROLS_INSIDE = ["PRPF31", "TMC4", "MBOAT7", "TTYH1"]
CONTROLS_OUTSIDE = ["PRKCG", "CACNG6", "PPP6R1"]

# Windows are reported only if they are at least this long. Shorter runs of
# unique sequence exist but do not survive a read length: a 150 bp read cannot
# sit inside an 80 bp window, so counting it as measurable depth would promise
# resolution the data does not have.
MIN_WINDOW = 300

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def read_fasta_region(reference: Path, chrom: str, start: int, end: int) -> str:
    """0-based half-open slice of one contig, via the .fai."""
    fai = {}
    for line in open(f"{reference}.fai"):
        f = line.split("\t")
        fai[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))
    if chrom not in fai:
        raise SystemExit(f"{chrom} not in {reference}.fai")
    length, offset, line_bases, line_width = fai[chrom]
    end = min(end, length)

    with open(reference, "rb") as fh:
        def seek_to(pos: int) -> int:
            return offset + (pos // line_bases) * line_width + (pos % line_bases)
        fh.seek(seek_to(start))
        raw = fh.read(seek_to(end) - seek_to(start))
    return raw.decode().replace("\n", "").upper()


def genes_from_gff(gff: Path, wanted: set[str]) -> dict[str, tuple[str, int, int, str]]:
    """``{symbol: (chrom, start0, end, strand)}`` for the gene features asked for."""
    found: dict[str, tuple[str, int, int, str]] = {}
    opener = gzip.open if str(gff).endswith(".gz") else open
    with opener(gff, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            attrs = dict(
                kv.split("=", 1) for kv in f[8].split(";") if "=" in kv
            )
            name = attrs.get("gene") or attrs.get("Name")
            if name not in wanted:
                continue
            start0, end = int(f[3]) - 1, int(f[4])
            # The union of records for one symbol, so a gene annotated in two
            # pieces is one interval rather than whichever piece came last.
            if name in found:
                c, s, e, st = found[name]
                if c == f[0]:
                    found[name] = (c, min(s, start0), max(e, end), st)
                    continue
            found[name] = (f[0], start0, end, f[6])
    return found


def ambiguous_positions(region: str, gene_offset: int, gene_len: int) -> list[bool]:
    """For each 100-mer start in the gene, whether a twin exists in the region.

    ``region`` is the whole LRC; ``gene_offset`` is where the gene begins inside
    it. A 100-mer is its own twin, so the position it came from is excluded --
    and so are overlapping positions, which differ from it by a shift rather than
    by paralogy.
    """
    # Every 25-mer of the region and its reverse complement, to its positions.
    seeds: dict[str, list[int]] = {}
    for i in range(len(region) - SEED + 1):
        seeds.setdefault(region[i:i + SEED], []).append(i)
    rc_region = revcomp(region)
    for i in range(len(rc_region) - SEED + 1):
        # Negative marks the reverse strand; the value is not a usable
        # coordinate, only a marker that a twin was found there.
        seeds.setdefault(rc_region[i:i + SEED], []).append(-(i + 1))

    out: list[bool] = []
    for g in range(gene_len - WINDOW + 1):
        start = gene_offset + g
        query = region[start:start + WINDOW]
        if "N" in query:
            out.append(True)        # unknowable, so not claimed as unique
            continue

        candidates: set[int] = set()
        for s in range(0, WINDOW - SEED + 1, SEED):
            for pos in seeds.get(query[s:s + SEED], ()):
                if pos >= 0:
                    cand = pos - s
                    if abs(cand - start) >= WINDOW:
                        candidates.add(cand)
                else:
                    candidates.add(pos)

        found_twin = False
        for cand in candidates:
            if cand >= 0:
                if cand < 0 or cand + WINDOW > len(region):
                    continue
                other = region[cand:cand + WINDOW]
            else:
                c = -cand - 1
                if c + WINDOW > len(rc_region):
                    continue
                other = rc_region[c:c + WINDOW]
                # A gene's own sequence appears on the reverse strand too; that
                # self-hit is a palindrome of itself, not a paralogue.
                rc_self = len(region) - (start + WINDOW)
                if abs(c - rc_self) < WINDOW:
                    continue
            mism = sum(1 for a, b in zip(query, other, strict=False) if a != b)
            if mism <= MAX_MISMATCH:
                found_twin = True
                break
        out.append(found_twin)
    return out


def unique_windows(ambiguous: list[bool], chrom: str, gene_start: int,
                   ) -> list[tuple[str, int, int]]:
    """Runs of unambiguous 100-mers, as merged intervals, longer than MIN_WINDOW.

    Merged, and that is not tidiness. A run of unambiguous 100-mer *starts*
    covers bases up to ``WINDOW`` past its last start, so two runs separated by a
    gap shorter than 100 bp produce intervals that overlap. Handed to
    `samtools depth` as two BED records, the overlap is counted twice, the
    gene's mean depth rises, and its copy number rises with it — the same class
    of bug as the doubled control loci that put λ₁ at 36.8 and halved every copy
    number in the cohort. Overlap is filtered here and asserted against in
    tests/test_loci_chm13.py.

    The union is the honest interval: every base in it is covered by at least
    one 100-mer with no twin in the LRC.
    """
    # Per base, not per run. A base is measurable if some unambiguous 100-mer
    # covers it; maximal runs of those bases are the windows. Defining it this
    # way makes overlap impossible by construction rather than by a merge, and —
    # the reason it was changed — it stops a merge from swallowing the ambiguous
    # stretches between two runs.
    #
    # That mattered at LILRA3, which is only 9.6% ambiguous and came out as one
    # window spanning the entire gene, the 9.6% included. Ambiguous positions
    # attract multi-mapping reads that fall below the MAPQ floor, so they
    # contribute near-zero depth to a MAPQ-20 mean: including them does not add
    # noise, it biases the gene's copy number *down* by roughly their fraction.
    n = len(ambiguous)
    covered = bytearray(n + WINDOW)
    for i, amb in enumerate(ambiguous):
        if not amb:
            for j in range(i, i + WINDOW):
                covered[j] = 1

    windows: list[tuple[str, int, int]] = []
    run_start = None
    for i, c in enumerate(bytes(covered) + b"\x00"):
        if c and run_start is None:
            run_start = i
        elif not c and run_start is not None:
            s, e = gene_start + run_start, gene_start + i
            if e - s >= MIN_WINDOW:
                windows.append((chrom, s, e))
            run_start = None
    return windows


def transfer_windows(reference: Path, grch38: Path, work: Path,
                     minimap2: str = "minimap2",
                     ) -> dict[str, list[tuple[str, int, int, float]]]:
    """Move `loci.UNIQUE_WINDOWS` onto CHM13 by aligning their sequence.

    Returns ``{gene: [(chrom, start0, end, identity), ...]}``. Identity is
    carried out so the caller can see how well each window transferred: these
    are two haplotypes of a polymorphic region, and a window that lands at 0.95
    rather than 0.99 is a different sequence being measured, not the same one at
    a new address.
    """
    import subprocess

    from lilrwgs import loci

    out: dict[str, list[tuple[str, int, int, float]]] = {}
    chrom, r_start, r_end = LRC_SEARCH
    region = read_fasta_region(reference, chrom, r_start, r_end)
    target = work / "chm13_lrc.fa"
    target.write_text(">lrc\n" + "\n".join(region[i:i + 60]
                                            for i in range(0, len(region), 60)) + "\n")

    for gene, windows in loci.UNIQUE_WINDOWS.items():
        placed: list[tuple[str, int, int, float]] = []
        for n, (wc, ws, we) in enumerate(windows, 1):
            seq = read_fasta_region(grch38, wc, ws, we)
            q = work / f"{gene}_{n}.fa"
            q.write_text(f">{gene}_{n}\n" +
                         "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")
            res = subprocess.run([minimap2, "-c", str(target), str(q)],
                                 capture_output=True, text=True)
            best = None
            for line in res.stdout.splitlines():
                f = line.split("\t")
                if len(f) < 12:
                    continue
                matches, block = int(f[9]), int(f[10])
                identity = matches / block if block else 0.0
                # Require most of the window to have transferred, so a short
                # high-identity hit to a paralogue cannot outrank the real one.
                if block < 0.8 * (we - ws):
                    continue
                if best is None or matches > best[0]:
                    best = (matches, int(f[7]), int(f[8]), identity)
            if best is None:
                print(f"  WARNING: {gene} window {n} ({we - ws} bp) did not "
                      "transfer to CHM13")
                continue
            _, t0, t1, identity = best
            placed.append((chrom, r_start + t0, r_start + t1, identity))
        out[gene] = placed
    return out


def locate_lilra3(reference: Path, lilra3_ref: Path,
                  minimap2: str = "minimap2") -> tuple[str, int, int, str, float]:
    """Find LILRA3 in CHM13 by aligning the calling reference to the LRC.

    Necessary because the annotation cannot help: `chm13v2.0_RefSeq_Liftoff` is
    lifted from the GRCh38 primary assembly, which carries the deletion and so
    has no LILRA3 to lift. Its absence there is a fact about the annotation's
    provenance, not about the sequence.

    The search space is the LILRB2-LILRA5 interval, which is where the gene sits
    in every LRC haplotype that has it, and which in CHM13 is 7,425 bp wider than
    in GRCh38 -- about one LILRA3.
    """
    import subprocess
    import tempfile

    chrom, r_start, r_end = LRC_SEARCH
    region = read_fasta_region(reference, chrom, r_start, r_end)

    with tempfile.TemporaryDirectory() as tmp:
        fa = Path(tmp) / "region.fa"
        fa.write_text(f">{chrom}_{r_start}\n" +
                      "\n".join(region[i:i + 60]
                                for i in range(0, len(region), 60)) + "\n")
        res = subprocess.run([minimap2, "-c", str(lilra3_ref), str(fa)],
                             capture_output=True, text=True)
    if res.returncode != 0:
        raise SystemExit(f"minimap2 failed: {res.stderr[:400]}")

    best = None
    for line in res.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        matches, block, mapq = int(f[9]), int(f[10]), int(f[11])
        identity = matches / block if block else 0.0
        # The real LILRA3 is a full-length, near-perfect hit; everything else in
        # this region is a partial paralogue at 59-87%. Requiring both length and
        # identity keeps a good LILRB2 hit from being mistaken for a poor LILRA3.
        if identity >= 0.98 and block >= 6_000 and mapq >= 30:
            if best is None or matches > best[0]:
                best = (matches, int(f[2]), int(f[3]), f[4], identity)

    if best is None:
        raise SystemExit(
            "no full-length LILRA3 match in the CHM13 LRC.\n"
            "  This reference does not carry the insertion allele, and the "
            "whole reason for using it does not hold. Stop rather than "
            "reporting LILRA3 from it.")
    _, q_start, q_end, strand, identity = best
    return chrom, r_start + q_start, r_start + q_end, strand, identity


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--reference", type=Path, required=True,
                   help="chm13v2.0.fa, indexed")
    p.add_argument("--grch38", type=Path,
                   default=Path("resources/reference/"
                                "GRCh38_full_analysis_set_plus_decoy_hla.fa"),
                   help="source of the validated unique-window sequence")
    p.add_argument("--gff", type=Path, required=True)
    p.add_argument("--lilra3-ref", type=Path,
                   default=Path("resources/bundle/references/LILRA3_named.fa"),
                   help="used only to report the LILRA3 interval for checking; "
                        "the interval itself is pinned in loci_chm13")
    p.add_argument("-o", "--output", type=Path, required=True)
    args = p.parse_args()

    wanted = set(LILR_GENES) | set(CONTROLS_INSIDE) | set(CONTROLS_OUTSIDE)
    genes = genes_from_gff(args.gff, wanted)
    missing = sorted(wanted - set(genes))
    if missing:
        print(f"warning: absent from the annotation: {', '.join(missing)}")

    chrom, r_start, r_end = LRC_SEARCH
    print(f"reading {chrom}:{r_start:,}-{r_end:,} ({r_end - r_start:,} bp) ...")
    region = read_fasta_region(args.reference, chrom, r_start, r_end)

    rows = []
    for name in sorted(genes):
        c, s, e, strand = genes[name]
        kind = ("control_inside" if name in CONTROLS_INSIDE else
                "control_outside" if name in CONTROLS_OUTSIDE else "gene")
        rows.append({"name": name, "kind": kind, "chrom": c,
                     "start": s, "end": e, "strand": strand, "note": ""})

    # LILRA3, which the annotation cannot supply. This is the whole reason for
    # the second reference, so it is derived and reported rather than assumed.
    a3c, a3s, a3e, a3strand, a3id = locate_lilra3(args.reference, args.lilra3_ref)
    print(f"LILRA3 located at {a3c}:{a3s + 1:,}-{a3e:,} "
          f"({a3e - a3s:,} bp, strand {a3strand}, identity {a3id:.4f})")
    rows.append({"name": "LILRA3", "kind": "gene", "chrom": a3c,
                 "start": a3s, "end": a3e, "strand": a3strand,
                 "note": f"located by alignment, identity {a3id:.4f}"})

    # The windows, transferred rather than re-derived. See the module docstring:
    # the documented criterion run against GRCh38 does not reproduce the windows
    # LILRA6's accuracy was validated on, so the validated sequence is moved
    # across and the criterion is used to audit the result.
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        placed = transfer_windows(args.reference, args.grch38, Path(tmp))

    # The audit: is each transferred window still free of twins in CHM13?
    ambiguity: dict[str, list[bool]] = {}
    for gene in ("LILRA6", "LILRB3"):
        if gene in genes:
            c, s, e, _ = genes[gene]
            print(f"scanning {gene} ({e - s:,} bp) for paralogue-unique "
                  "100-mers, as an audit ...")
            ambiguity[gene] = ambiguous_positions(region, s - r_start, e - s)
            print(f"  {sum(ambiguity[gene]) / max(len(ambiguity[gene]), 1):.1%} "
                  "of 100-mers have a twin")

    for gene, windows in placed.items():
        strand = genes[gene][3] if gene in genes else "."
        for i, (wc, ws, we, identity) in enumerate(windows, 1):
            flag = ""
            if gene in ambiguity and gene in genes:
                g_start = genes[gene][1]
                lo, hi = ws - g_start, we - g_start - WINDOW + 1
                inside = ambiguity[gene][max(0, lo):max(0, hi)]
                if inside:
                    amb_frac = sum(inside) / len(inside)
                    flag = f", {amb_frac:.1%} ambiguous by the 100-mer audit"
                    if amb_frac > 0.10:
                        print(f"  WARNING: {gene} window {i} is {amb_frac:.1%} "
                              "ambiguous in CHM13 — it is not measuring what it "
                              "measured on GRCh38")
            rows.append({"name": f"{gene}_unique_{i}", "kind": "unique_window",
                         "chrom": wc, "start": ws, "end": we, "strand": strand,
                         "note": f"{we - ws} bp, transferred at identity "
                                 f"{identity:.4f}{flag}"})

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as fh:
        fh.write("name\tkind\tchrom\tstart\tend\tstrand\tnote\n")
        for r in rows:
            fh.write(f"{r['name']}\t{r['kind']}\t{r['chrom']}\t{r['start']}\t"
                     f"{r['end']}\t{r['strand']}\t{r['note']}\n")
    print(f"\nwrote {len(rows)} intervals -> {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
