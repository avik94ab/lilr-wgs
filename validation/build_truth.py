#!/usr/bin/env python3
"""Build a truth set from the HPRC panels that the pipeline aligns against.

The panel FASTA headers carry their donor and haplotype:

    >HG00097|hap1|LILRA6|copy1|HPRC|EUR|len=6358

so the same files are both the alignment target and a truth set. For a donor
that also has a 1000 Genomes 30x CRAM, the assembly gives two answers the
pipeline can be scored on:

* **copy number** — count `copyN` entries per donor per gene across both
  haplotypes;
* **allele sequence** — the assembly contig itself, base for base.

231 HPRC donors are represented and 101 of them are in the 1000 Genomes 2504
panel, which is the validation set.

**The circularity, stated plainly.** A donor in that overlap is being aligned
against a panel containing its own haplotypes, which inflates recruitment and
copy-number accuracy relative to an unseen sample. `--leave-one-donor-out` writes
panels with a donor's own entries removed, and the score from those is the one
that generalises. Reporting only the as-is number would be reporting the easy
case.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def parse_header(header: str) -> dict | None:
    """`HG00097|hap1|LILRA6|copy1|HPRC|EUR|len=6358` -> its fields.

    Returns None for a header that does not match, rather than guessing: a panel
    built to a different convention should be noticed, not silently
    half-parsed.
    """
    parts = header.lstrip(">").strip().split("|")
    if len(parts) < 4:
        return None
    donor, hap, gene, copy = parts[0], parts[1], parts[2], parts[3]
    if not hap.startswith("hap") or not copy.startswith("copy"):
        return None
    return {"donor": donor, "hap": hap, "gene": gene, "copy": copy,
            "population": parts[5] if len(parts) > 5 else ""}


def read_panels(panel_dir: Path) -> dict[str, list[dict]]:
    """Every panel entry, keyed by gene, with its sequence."""
    by_gene: dict[str, list[dict]] = defaultdict(list)
    for fasta in sorted(panel_dir.glob("*.fasta")):
        gene = fasta.stem
        entry, seq = None, []
        for line in fasta.read_text().splitlines():
            if line.startswith(">"):
                if entry:
                    entry["seq"] = "".join(seq)
                    by_gene[gene].append(entry)
                entry, seq = parse_header(line), []
            elif entry:
                seq.append(line.strip())
        if entry:
            entry["seq"] = "".join(seq)
            by_gene[gene].append(entry)
    return by_gene


# A donor present in at least this many of the 11 panels is taken to have a
# well-resolved LRC, so absence from one more panel means the gene is deleted
# rather than unassembled.
MIN_GENES_FOR_ABSENCE = 8


def copy_number_truth(by_gene: dict[str, list[dict]],
                      ) -> tuple[dict[str, dict[str, int]], dict[str, set[str]]]:
    """Per donor, per gene, the number of copies across both haplotypes.

    Absence from a gene's panel is ambiguous in general — it could be a true
    deletion or an assembly that did not resolve the locus — and calling it copy
    number 0 unconditionally would manufacture deletions, which is exactly the
    error the pipeline itself is careful to avoid.

    But leaving it out entirely loses the most interesting class in the dataset.
    LILRA3's panel holds 352 sequences from 208 donors, against 232 donors
    represented overall: the missing 24 are deletion homozygotes, and dropping
    them would leave the truth set with no CN 0 at the one gene where CN 0 is
    common. So absence is read as zero *only* for a donor whose LRC is otherwise
    well resolved, and which entries were inferred rather than counted is
    returned alongside so a consumer can score with or without them.
    """
    truth: dict[str, dict[str, int]] = defaultdict(dict)
    for gene, entries in by_gene.items():
        counts: dict[str, int] = defaultdict(int)
        for e in entries:
            counts[e["donor"]] += 1
        for donor, n in counts.items():
            truth[donor][gene] = n

    inferred: dict[str, set[str]] = defaultdict(set)
    all_genes = set(by_gene)
    for donor, genes in truth.items():
        if len(genes) < MIN_GENES_FOR_ABSENCE:
            continue
        for gene in all_genes - set(genes):
            truth[donor][gene] = 0
            inferred[donor].add(gene)
    return truth, inferred


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--panels", type=Path, default=Path("resources/gdna"))
    p.add_argument("--kgp-samples", type=Path,
                   help="file of 1000 Genomes sample IDs, to compute the overlap")
    p.add_argument("--outdir", type=Path, default=Path("validation/truth"))
    p.add_argument("--leave-one-donor-out", type=Path,
                   help="write donor-excluded panels under this directory")
    args = p.parse_args()

    by_gene = read_panels(args.panels)
    if not by_gene:
        raise SystemExit(f"no panels found under {args.panels}")
    truth, inferred = copy_number_truth(by_gene)
    args.outdir.mkdir(parents=True, exist_ok=True)

    genes = sorted(by_gene)
    with (args.outdir / "copy_number.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["donor", "gene", "copies", "source"])
        for donor in sorted(truth):
            for gene in genes:
                if gene in truth[donor]:
                    source = ("inferred_absent" if gene in inferred.get(donor, ())
                              else "counted")
                    w.writerow([donor, gene, truth[donor][gene], source])

    overlap = []
    if args.kgp_samples and args.kgp_samples.exists():
        kgp = {ln.strip() for ln in args.kgp_samples.read_text().splitlines() if ln.strip()}
        overlap = sorted(set(truth) & kgp)
        (args.outdir / "overlap_1kgp.txt").write_text("\n".join(overlap) + "\n")

    if args.leave_one_donor_out and overlap:
        for donor in overlap:
            donor_dir = args.leave_one_donor_out / donor
            donor_dir.mkdir(parents=True, exist_ok=True)
            for gene, entries in by_gene.items():
                kept = [e for e in entries if e["donor"] != donor]
                (donor_dir / f"{gene}.fasta").write_text("".join(
                    f">{e['donor']}|{e['hap']}|{e['gene']}|{e['copy']}\n{e['seq']}\n"
                    for e in kept))

    summary = {
        "n_donors": len(truth),
        "n_genes": len(genes),
        "n_overlap_1kgp": len(overlap),
        "n_inferred_absent": sum(len(v) for v in inferred.values()),
        "entries_per_gene": {g: len(e) for g, e in sorted(by_gene.items())},
        "cn_distribution": {
            g: dict(sorted(
                (n, sum(1 for d in truth.values() if d.get(g) == n))
                for n in sorted({d[g] for d in truth.values() if g in d})))
            for g in genes
        },
    }
    (args.outdir / "truth_summary.json").write_text(json.dumps(summary, indent=2))

    print(f"{len(truth)} donors, {len(genes)} genes -> {args.outdir}")
    if overlap:
        print(f"{len(overlap)} donors also have a 1000 Genomes 30x CRAM")
    for gene in genes:
        print(f"  {gene}: {summary['cn_distribution'][gene]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
