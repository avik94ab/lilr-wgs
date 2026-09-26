#!/usr/bin/env python3
"""Score called protein, coding sequence and haplotypes against HPRC truth.

`docs/variant_calling.md` §2.4: the step that licenses a claim about allele
sequence. `compare_cn.py` scores copy number; this scores what the genotype stage
emits, with the protein first because that is the claim a user of the output
makes.

Truth is the assembly itself. Each donor's panel entries (`resources/gdna`, the
same files recruitment aligns against) are annotated exactly as the pipeline
annotates a called haplotype -- miniprot against the gene's reference protein,
then :func:`lilrwgs.sequences.extract_sequences` -- so a disagreement is between
sequences, not between two annotation methods.

Three levels, each scored per gene and never as one aggregate:

* **protein** -- each called haplotype's protein against the truth haplotype it
  matches best (one-to-one). `exact`, `consistent` (every resolved residue
  right, some X) or `wrong`. Amino-acid changes against the reference protein
  are scored as recovered / masked / missed / wrong, and called changes as
  correct or false. A *reference baseline* is printed beside every rate: the
  share of truth haplotypes whose protein *is* the reference, which is what
  emitting the reference for everyone would score. An exact rate means nothing
  without it.
* **coding sequence** -- the same, base by base on the cDNA, so synonymous
  changes count.
* **gDNA** -- truth aligned to the calling reference; SNV genotype concordance
  at callable positions and haplotype base errors by region (CDS, intron,
  flank, and for LILRB1/LILRB4 the shared block), plus switch errors between
  consecutive heterozygous sites.

The called gDNA is SNV-only by construction (indels never reach the consensus),
so a truth indel is not scorable as a base and is counted separately; in the
protein it is scored as it lands.

Run it twice -- once on the as-is run and once on the leave-one-donor-out run --
and report both: a donor in the overlap is recruited against a panel containing
its own haplotypes, and allele sequence runs through recruitment.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import os
import pickle
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "validation"))

from build_truth import parse_header                          # noqa: E402
from lilrwgs.genotype import _parse_gff                        # noqa: E402
from lilrwgs.sequences import extract_sequences, stop_codon_status  # noqa: E402

DEFAULT_GENES = ["LILRA1", "LILRA2", "LILRA3", "LILRA4", "LILRA5",
                 "LILRB1", "LILRB2", "LILRB4", "LILRB5"]
SHARED_GENES = {"LILRB1", "LILRB4"}       # config shared_groups, minus LILRA6/LILRB3
_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


# --------------------------------------------------------------------------
# inputs
# --------------------------------------------------------------------------

def read_fasta(path: Path) -> list[tuple[str, str]]:
    out, hdr, cur = [], None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if hdr is not None:
                out.append((hdr, "".join(cur)))
            hdr, cur = line[1:].strip(), []
        else:
            cur.append(line.strip())
    if hdr is not None:
        out.append((hdr, "".join(cur)))
    return out


def one_seq(path: Path) -> str:
    recs = read_fasta(path) if path.exists() else []
    return recs[0][1].upper() if recs else ""


def revcomp(s: str) -> str:
    return s.translate(_RC)[::-1]


def truth_haplotypes(panel: Path, donors: set[str]) -> dict[str, list[dict]]:
    """Every panel entry for these donors, in header order."""
    out: dict[str, list[dict]] = defaultdict(list)
    for hdr, seq in read_fasta(panel):
        h = parse_header(hdr)
        if h and h["donor"] in donors:
            out[h["donor"]].append({"label": f"{h['hap']}.{h['copy']}",
                                    "seq": seq.upper()})
    return out


# --------------------------------------------------------------------------
# annotation and alignment (the slow parts; cached across runs)
# --------------------------------------------------------------------------

def annotate(seq: str, protein_faa: str, miniprot: str = "miniprot") -> dict:
    """Protein, cDNA and stop status, the way genotype._finish_haplotype does."""
    with tempfile.TemporaryDirectory() as tmp:
        fa = Path(tmp) / "h.fa"
        fa.write_text(f">h\n{seq}\n")
        res = subprocess.run([miniprot, "--gff", str(fa), protein_faa],
                             capture_output=True, text=True, timeout=300)
    coords = _parse_gff(res.stdout)
    if not coords:
        return {"cdna": "", "protein": "", "stop": "", "coords": None}
    _, cdna, protein = extract_sequences(seq, coords)
    return {"cdna": cdna, "protein": protein,
            "stop": stop_codon_status(cdna, coords), "coords": coords}


def orient(seq: str, ref: str, k: int = 25) -> str:
    """The truth haplotype on the reference's strand. Panels for minus-strand
    genes are stored in gene orientation; the calling references are not."""
    kmers = {ref[i:i + k] for i in range(0, len(ref) - k, 50)}
    fwd = sum(seq[i:i + k] in kmers for i in range(0, len(seq) - k))
    rc = revcomp(seq)
    rev = sum(rc[i:i + k] in kmers for i in range(0, len(rc) - k))
    return seq if fwd >= rev else rc


def align_to_ref(ref: str, seq: str) -> dict:
    """Map a truth haplotype onto reference coordinates.

    Returns ``base``: one character per reference position (the truth base, or
    '-' where truth deletes it, or '' where the truth sequence does not reach),
    and ``ins``: reference positions after which truth inserts bases.
    """
    from Bio import Align

    a = Align.PairwiseAligner()
    a.mode = "global"
    a.match_score, a.mismatch_score = 2, -3
    a.open_gap_score, a.extend_gap_score = -7, -1
    a.end_gap_score = 0
    aln = a.align(ref, seq)[0]
    base = [""] * len(ref)
    ins: set[int] = set()
    blocks = list(zip(aln.aligned[0], aln.aligned[1]))
    for (rs, re_), (qs, qe) in blocks:
        for i in range(re_ - rs):
            base[rs + i] = seq[qs + i]
    for ((_, r_end), (_, q_end)), ((r_next, _), (q_next, _)) in zip(blocks, blocks[1:]):
        if r_next > r_end:                    # truth lacks reference bases
            for p in range(r_end, r_next):
                base[p] = "-"
        if q_next > q_end:                    # truth has extra bases
            ins.add(r_end - 1)
    return {"base": base, "ins": sorted(ins)}


def _truth_job(args):
    gene, donor, label, seq, ref, faa, miniprot = args
    ann = annotate(seq, faa, miniprot)
    ann.pop("coords")
    aligned = align_to_ref(ref, orient(seq, ref))
    return (gene, donor, label), {**ann, **aligned}


def build_truth(genes, donors, panels: Path, refs: Path, exon_dir: Path,
                cache_path: Path | None, jobs: int, miniprot: str) -> dict:
    cache = {}
    if cache_path and cache_path.exists():
        cache = pickle.loads(cache_path.read_bytes())
    todo, keys = [], {}
    for gene in genes:
        ref = one_seq(refs / f"{gene}_named.fa")
        faa = str(exon_dir / f"{gene}_protein.faa")
        for donor, haps in truth_haplotypes(panels / f"{gene}.fasta", donors).items():
            for h in haps:
                digest = hashlib.sha1(f"{gene}|{h['seq']}".encode()).hexdigest()
                keys[(gene, donor, h["label"])] = digest
                if digest not in cache:
                    todo.append((digest, (gene, donor, h["label"], h["seq"], ref,
                                          faa, miniprot)))
    if todo:
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            for (digest, _), (_, res) in zip(todo, ex.map(_truth_job,
                                                           [t for _, t in todo])):
                cache[digest] = res
        if cache_path:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            cache_path.write_bytes(pickle.dumps(cache))
    truth: dict = defaultdict(lambda: defaultdict(list))
    for (gene, donor, label), digest in sorted(keys.items()):
        truth[gene][donor].append({"label": label, **cache[digest]})
    return truth


# --------------------------------------------------------------------------
# comparison
# --------------------------------------------------------------------------

def residue_errors(called: str, truth: str) -> int:
    """Resolved residues that disagree. X is unknown, never wrong. A length
    difference counts every unmatched residue: the pipeline cannot call indels,
    so a truth indel is a real miss in the protein it reports."""
    if len(called) == len(truth):
        return sum(1 for c, t in zip(called, truth) if c != "X" and c != t)
    from Bio import Align
    a = Align.PairwiseAligner()
    a.mode = "global"
    a.match_score, a.mismatch_score = 1, -1
    a.open_gap_score, a.extend_gap_score = -3, -1
    aln = a.align(called, truth)[0]
    errs, n_aligned_called = 0, 0
    for (cs, ce), (ts, te) in zip(*aln.aligned):
        for i in range(ce - cs):
            c, t = called[cs + i], truth[ts + i]
            errs += c != "X" and c != t
        n_aligned_called += ce - cs
    unmatched = (len(called) - n_aligned_called) + (len(truth) - n_aligned_called)
    return errs + unmatched


def classify_protein(called: str, truth: str) -> str:
    if not called:
        return "missing"
    if called == truth:
        return "exact"
    if len(called) == len(truth) and residue_errors(called, truth) == 0:
        return "consistent"
    return "wrong"


def best_assignment(called: list, truth: list, cost) -> tuple[tuple[int, ...], int]:
    """One-to-one matching of called to truth haplotypes minimising ``cost``."""
    best, best_cost = None, None
    for perm in itertools.permutations(range(len(truth)), len(called)):
        c = sum(cost(called[i], truth[j]) for i, j in enumerate(perm))
        if best_cost is None or c < best_cost:
            best, best_cost = perm, c
    return best, best_cost


def variant_tally(called: str, truth: str, ref: str, masked: str = "X") -> Counter:
    """Changes against the reference, at one aligned haplotype pair.

    ``masked`` is the unknown symbol: X in a protein, N in DNA (where N in a
    protein is asparagine, so the two cannot share a default). Only defined
    where all three have the same length; a length change is reported by the
    protein class instead.
    """
    t = Counter()
    if not (len(called) == len(truth) == len(ref)):
        t["unscorable_length"] += 1
        return t
    for c, tr, r in zip(called, truth, ref):
        if tr != r:
            t["truth_changes"] += 1
            if c == tr:
                t["recovered"] += 1
            elif c == masked:
                t["masked"] += 1
            elif c == r:
                t["missed_ref"] += 1
            else:
                t["wrong_alt"] += 1
        if c not in (masked, r):
            t["called_changes"] += 1
            t["called_correct" if c == tr else "called_false"] += 1
        t["positions"] += 1
        t["x"] += c == masked
    return t


def unordered_match(called: list[str], truth: list[str]) -> str:
    """Per-position multisets, phase ignored: is every resolved position's set
    of residues right? Separates a phase error from a residue error."""
    if len({len(s) for s in called + truth}) != 1:
        return "length"
    for col in zip(*called, *truth):
        c, t = col[:len(called)], col[len(called):]
        if "X" in c:
            continue
        if sorted(c) != sorted(t):
            return "residue"
    return "ok"


# --------------------------------------------------------------------------
# gDNA
# --------------------------------------------------------------------------

def region_classes(gff: Path, length: int) -> list[str]:
    coords = _parse_gff(gff.read_text())
    cls = ["flank"] * length
    if coords:
        for p in range(coords["mrna_start"], coords["mrna_end"]):
            cls[p] = "intron"
        for s, e in coords["exons"]:
            for p in range(s, e):
                cls[p] = "cds"
    return cls


def shared_positions(track: Path) -> set[int]:
    if not track.exists():
        return set()
    out = set()
    with track.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if float(row["shared_fraction"]) > 0.5:
                out.add(int(row["pos"]) - 1)
    return out


def score_gdna(called: list[str], truth: list[dict], ref: str, cls: list[str],
               shared: set[int]) -> dict:
    """Genotype and haplotype concordance at single-base positions."""
    k = len(called)
    tb = [t["base"] for t in truth]
    ins = set().union(*[set(t["ins"]) for t in truth]) if truth else set()
    out: dict = defaultdict(Counter)

    def region(p):
        return "shared" if p in shared else cls[p]

    # Haplotype assignment on gDNA mismatches, independent of the protein's.
    def hap_cost(c, t):
        return sum(1 for p in range(len(ref))
                   if c[p] not in ("N",) and t[p] not in ("", "-") and c[p] != t[p])
    perm, _ = best_assignment(called, tb, hap_cost)

    het_orientation = []
    for p in range(len(ref)):
        tcol = [t[p] for t in tb]
        if "" in tcol:
            continue
        r = region(p)
        if "-" in tcol or p in ins:
            out[r]["truth_indel_positions"] += 1
            continue
        out[r]["positions"] += 1
        ccol = [c[p] for c in called]
        truth_var = any(b != ref[p] for b in tcol)
        if "N" in ccol:
            out[r]["uncallable"] += 1
            if truth_var:
                out[r]["truth_var_masked"] += 1
            continue
        out[r]["callable"] += 1
        gt_ok = sorted(ccol) == sorted(tcol)
        out[r]["gt_ok"] += gt_ok
        called_var = any(b != ref[p] for b in ccol)
        if truth_var:
            out[r]["truth_var"] += 1
            out[r]["truth_var_gt_ok"] += gt_ok
        if called_var:
            out[r]["called_var"] += 1
            out[r]["called_var_gt_ok"] += gt_ok
        for i, j in enumerate(perm):
            out[r]["hap_bases"] += 1
            out[r]["hap_errors"] += called[i][p] != tb[j][p]
        if k == 2 and tcol[0] != tcol[1] and gt_ok:
            het_orientation.append(called[perm.index(0)][p] == tcol[0])
    switches = sum(1 for a, b in zip(het_orientation, het_orientation[1:]) if a != b)
    out["all"]["het_sites"] = len(het_orientation)
    out["all"]["switches"] = switches
    return out


# --------------------------------------------------------------------------
# per-run scoring
# --------------------------------------------------------------------------

def score_run(run: Path, genes, truth, refs: Path, exon_dir: Path,
              samples: list[str]) -> tuple[list[dict], dict]:
    rows, gdna = [], defaultdict(lambda: defaultdict(Counter))
    for gene in genes:
        ref = one_seq(refs / f"{gene}_named.fa")
        ref_protein = one_seq(exon_dir / f"{gene}_protein.faa")
        ref_ann = extract_sequences(ref, _parse_gff((exon_dir / f"{gene}.gff").read_text()))
        ref_cdna = ref_ann[1]
        cls = region_classes(exon_dir / f"{gene}.gff", len(ref))
        for donor in samples:
            thaps = truth[gene].get(donor, [])
            marker = run / "genotypes" / "markers" / f"{donor}__{gene}.json"
            m = json.loads(marker.read_text()) if marker.exists() else {}
            cn = m.get("cn")
            gdir = run / "genotypes" / donor
            haps = [h["hap"] for h in m.get("haplotypes", [])]
            called = [{"hap": h,
                       "protein": one_seq(gdir / f"{gene}_hap{h}_protein.fasta"),
                       "cdna": one_seq(gdir / f"{gene}_hap{h}_cdna.fasta"),
                       "gdna": one_seq(gdir / f"{gene}_hap{h}_gdna.fasta")}
                      for h in haps]
            base = {"gene": gene, "donor": donor, "cn": cn,
                    "truth_cn": len(thaps), "status": m.get("status", "missing")}
            if len(called) != len(thaps) or not called:
                rows.append({**base, "hap": "", "class": "cn_mismatch"
                             if len(called) != len(thaps) else "no_haplotypes"})
                continue
            tprot = [t["protein"] for t in thaps]
            perm, _ = best_assignment(
                [c["protein"] for c in called], tprot,
                lambda c, t: residue_errors(c, t) if c else 10**6)
            unordered = unordered_match([c["protein"] for c in called],
                                        [tprot[j] for j in perm])
            for i, j in enumerate(perm):
                c, t = called[i], thaps[j]
                vt = variant_tally(c["protein"], t["protein"], ref_protein)
                ct = variant_tally(c["cdna"], t["cdna"], ref_cdna, masked="N")
                rows.append({
                    **base, "hap": c["hap"], "truth_hap": t["label"],
                    "class": classify_protein(c["protein"], t["protein"]),
                    "residue_errors": residue_errors(c["protein"], t["protein"])
                    if c["protein"] else "",
                    "protein_len": len(c["protein"]),
                    "truth_len": len(t["protein"]),
                    "x": c["protein"].count("X"),
                    "truth_is_ref": t["protein"] == ref_protein,
                    "truth_stop": t["stop"],
                    "unordered": unordered,
                    **{f"aa_{k}": v for k, v in vt.items()},
                    **{f"nt_{k}": v for k, v in ct.items()},
                    "diffs": ";".join(
                        f"{p + 1}{ref_protein[p] if p < len(ref_protein) else '?'}"
                        f">{t['protein'][p]}/{c['protein'][p]}"
                        for p in range(min(len(c["protein"]), len(t["protein"])))
                        if c["protein"][p] not in ("X", t["protein"][p]))
                    if len(c["protein"]) == len(t["protein"]) else "length",
                })
            if all(c["gdna"] and len(c["gdna"]) == len(ref) for c in called):
                g = score_gdna([c["gdna"] for c in called], thaps, ref, cls,
                               shared_positions(gdir / f"{gene}.callability.tsv")
                               if gene in SHARED_GENES else set())
                key = (gene, cn)
                for region, cnt in g.items():
                    gdna[key][region].update(cnt)
    return rows, gdna


# --------------------------------------------------------------------------
# report
# --------------------------------------------------------------------------

def pct(a, b):
    return f"{100 * a / b:5.1f}%" if b else "    - "


def report(rows, gdna, run_label: str, genes, n_donors: int) -> str:
    L = []
    w = L.append
    w(f"# Sequence accuracy against HPRC assemblies: {run_label}")
    w("")
    w(f"{n_donors} donors x {len(genes)} genes ({', '.join(genes)}).")
    w("Truth: each donor's own assembly haplotypes, annotated with the pipeline's")
    w("own miniprot + extract_sequences. Diagnostic of the pipeline as run; the")
    w("in-panel run is optimistic (see docs/variant_calling.md s2.4).")
    w("")
    hap_rows = [r for r in rows if r.get("hap") != ""]
    other = [r for r in rows if r.get("hap") == ""]
    if other:
        w(f"NOT SCORED ({len(other)}): " + ", ".join(
            f"{r['donor']}/{r['gene']} {r['class']} (cn {r['cn']} vs truth "
            f"{r['truth_cn']})" for r in other))
        w("")

    def groups():
        for gene in genes:
            gr = [r for r in hap_rows if r["gene"] == gene]
            cns = sorted({r["cn"] for r in gr})
            if len(cns) > 1:
                for cn in cns:
                    yield f"{gene} CN{cn}", [r for r in gr if r["cn"] == cn]
            yield gene if len(cns) <= 1 else f"{gene} all", gr

    w("## 1. Proteins, per haplotype (called vs the truth haplotype it matches)")
    w("")
    w("exact = identical; consistent = every resolved residue right, some X;")
    w("wrong = at least one resolved residue wrong, or a length difference.")
    w("ref baseline = truth haplotypes whose protein IS the reference protein:")
    w("what emitting the reference for everyone would score as exact.")
    w("")
    w(f"{'gene':<12}{'haps':>5}{'exact':>9}{'consist.':>9}{'correct*':>9}"
      f"{'wrong':>8}{'ref base':>9}{'X res':>7}")
    tot = Counter()
    for name, gr in groups():
        c = Counter(r["class"] for r in gr)
        n = len(gr)
        xs = sum(r["x"] for r in gr)
        resid = sum(r["protein_len"] for r in gr)
        base = sum(bool(r["truth_is_ref"]) for r in gr)
        w(f"{name:<12}{n:>5}{pct(c['exact'], n):>9}{pct(c['consistent'], n):>9}"
          f"{pct(c['exact'] + c['consistent'], n):>9}{pct(c['wrong'] + c['missing'], n):>8}"
          f"{pct(base, n):>9}{pct(xs, resid):>7}")
        if " CN" not in name:
            tot.update({"n": n, "exact": c["exact"], "consistent": c["consistent"],
                        "wrong": c["wrong"] + c["missing"], "base": base,
                        "x": xs, "res": resid})
    w(f"{'all 9':<12}{tot['n']:>5}{pct(tot['exact'], tot['n']):>9}"
      f"{pct(tot['consistent'], tot['n']):>9}"
      f"{pct(tot['exact'] + tot['consistent'], tot['n']):>9}"
      f"{pct(tot['wrong'], tot['n']):>8}{pct(tot['base'], tot['n']):>9}"
      f"{pct(tot['x'], tot['res']):>7}")
    w("(* correct = exact + consistent)")
    w("")

    w("## 2. Amino-acid changes against the reference protein")
    w("")
    w("Of the changes truth carries: recovered / masked (X) / missed (called")
    w("reference) / wrong residue. Of the changes called: correct / false.")
    w("")
    w(f"{'gene':<12}{'truth':>7}{'recov.':>8}{'masked':>8}{'missed':>8}{'wrong':>7}"
      f"{'called':>8}{'correct':>9}{'false':>7}")
    for name, gr in groups():
        s = Counter()
        for r in gr:
            s.update({k[3:]: v for k, v in r.items()
                      if k.startswith("aa_") and isinstance(v, int)})
        tc, cc = s["truth_changes"], s["called_changes"]
        w(f"{name:<12}{tc:>7}{pct(s['recovered'], tc):>8}{pct(s['masked'], tc):>8}"
          f"{pct(s['missed_ref'], tc):>8}{pct(s['wrong_alt'], tc):>7}{cc:>8}"
          f"{pct(s['called_correct'], cc):>9}{pct(s['called_false'], cc):>7}"
          + (f"  ({s['unscorable_length']} hap(s) length-changed)"
             if s["unscorable_length"] else ""))
    w("")

    w("## 3. Protein genotype per donor (phase ignored)")
    w("")
    w("ok = every resolved position has the right residues on the right number")
    w("of haplotypes. A donor with wrong proteins but an 'ok' genotype has a")
    w("phase error, not a residue error.")
    w("")
    w(f"{'gene':<12}{'donors':>7}{'all haps right':>16}{'phase only':>12}"
      f"{'residue':>9}{'length':>8}")
    for name, gr in groups():
        by = defaultdict(list)
        for r in gr:
            by[r["donor"]].append(r)
        n = len(by)
        right = sum(all(r["class"] in ("exact", "consistent") for r in v)
                    for v in by.values())
        phase = sum(not all(r["class"] in ("exact", "consistent") for r in v)
                    and v[0]["unordered"] == "ok" for v in by.values())
        resid = sum(v[0]["unordered"] == "residue" for v in by.values())
        length = sum(v[0]["unordered"] == "length" for v in by.values())
        w(f"{name:<12}{n:>7}{pct(right, n):>16}{phase:>12}{resid:>9}{length:>8}")
    w("")

    w("## 4. Coding sequence (cDNA), per haplotype, changes against the reference")
    w("")
    w(f"{'gene':<12}{'bases':>8}{'X/N':>7}{'truth':>7}{'recov.':>8}{'masked':>8}"
      f"{'missed':>8}{'called':>8}{'false':>7}")
    for name, gr in groups():
        s = Counter()
        for r in gr:
            s.update({k[3:]: v for k, v in r.items()
                      if k.startswith("nt_") and isinstance(v, int)})
        tc, cc = s["truth_changes"], s["called_changes"]
        w(f"{name:<12}{s['positions']:>8}{pct(s['x'], s['positions']):>7}{tc:>7}"
          f"{pct(s['recovered'], tc):>8}{pct(s['masked'], tc):>8}"
          f"{pct(s['missed_ref'] + s['wrong_alt'], tc):>8}{cc:>8}"
          f"{pct(s['called_false'], cc):>7}")
    w("")
    w("(cDNA N counted in X/N; missed = called reference or a wrong base)")
    w("")

    w("## 5. gDNA by region: genotype and haplotype concordance")
    w("")
    w("callable = share of scorable positions called; GT@var = genotype right at")
    w("truth-variant positions that were called; false = called-variant")
    w("positions whose genotype is wrong; hap err = haplotype base errors per")
    w("10 kb called; indel = positions under a truth indel (unscorable).")
    w("")
    w(f"{'gene':<12}{'region':<8}{'callable':>9}{'GT@var':>9}{'var N':>7}"
      f"{'false':>8}{'hap err/10kb':>13}{'indel':>7}")
    for (gene, cn), regs in sorted(gdna.items(), key=lambda kv: (kv[0][0], kv[0][1] or 0)):
        label = gene + (f" CN{cn}" if gene == "LILRA3" else "")
        for region in ("cds", "intron", "flank", "shared"):
            s = regs.get(region)
            if not s or not s["positions"]:
                continue
            tv = s["truth_var"] + s["truth_var_masked"]
            w(f"{label:<12}{region:<8}{pct(s['callable'], s['positions']):>9}"
              f"{pct(s['truth_var_gt_ok'], s['truth_var']):>9}{pct(s['truth_var_masked'], tv):>7}"
              f"{pct(s['called_var'] - s['called_var_gt_ok'], s['called_var']):>8}"
              f"{(10_000 * s['hap_errors'] / s['hap_bases']) if s['hap_bases'] else 0:>13.1f}"
              f"{s['truth_indel_positions']:>7}")
        a = regs.get("all", Counter())
        if a.get("het_sites"):
            w(f"{label:<12}{'phasing':<8} {a['switches']} switches over "
              f"{a['het_sites']} correctly genotyped het sites")
    w("")

    wrong = [r for r in hap_rows if r["class"] in ("wrong", "missing")]
    w(f"## 6. Every wrong protein ({len(wrong)})")
    w("")
    w("diffs: position ref>truth/called, resolved residues only")
    for r in sorted(wrong, key=lambda r: (r["gene"], r["donor"], r["hap"])):
        w(f"{r['gene']:<8}{r['donor']:<9}hap{r['hap']} vs {r['truth_hap']:<12}"
          f"errors={r['residue_errors']:<3} len {r['protein_len']}/{r['truth_len']} "
          f"genotype={r['unordered']:<8} {r['diffs']}")
    return "\n".join(L) + "\n"


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("run", type=Path, help="pipeline outdir (has genotypes/)")
    p.add_argument("--label", default=None)
    p.add_argument("--panels", type=Path, default=REPO / "resources" / "gdna",
                   help="TRUTH panels -- always the full ones, never LODO")
    p.add_argument("--truth-cn", type=Path,
                   default=REPO / "validation" / "truth" / "copy_number.tsv")
    p.add_argument("--refs", type=Path,
                   default=REPO / "resources" / "bundle" / "references")
    p.add_argument("--exon-coords", type=Path,
                   default=REPO / "resources" / "bundle" / "exon_coords")
    p.add_argument("--genes", nargs="+", default=DEFAULT_GENES)
    p.add_argument("--cache", type=Path, default=None,
                   help="pickle of truth annotation/alignment, reused across runs")
    p.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    p.add_argument("--miniprot", default="miniprot")
    p.add_argument("-o", "--out", type=Path, required=True)
    p.add_argument("--tsv", type=Path, help="per-haplotype table")
    args = p.parse_args(argv)

    truth_donors = set()
    with args.truth_cn.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            truth_donors.add(row["donor"])
    called = {m.name.split("__")[0]
              for m in (args.run / "genotypes" / "markers").glob("*.json")}
    samples = sorted(truth_donors & called)
    truth = build_truth(args.genes, set(samples), args.panels, args.refs,
                        args.exon_coords, args.cache, args.jobs, args.miniprot)
    rows, gdna = score_run(args.run, args.genes, truth, args.refs,
                           args.exon_coords, samples)
    text = report(rows, gdna, args.label or str(args.run), args.genes, len(samples))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(text)
    if args.tsv:
        cols = sorted({k for r in rows for k in r},
                      key=lambda k: (k not in ("gene", "donor", "hap"), k))
        with args.tsv.open("w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            wr.writeheader()
            wr.writerows(rows)
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
