#!/usr/bin/env python3
"""Build a lilr-wgs manifest from a 1000 Genomes sequence index.

The published index is a 22-column TSV whose first column is an ``ftp://`` path
at ENA. EBI serves the same tree over HTTPS, and htslib can read a CRAM there
directly — so the manifest carries HTTPS URLs and nothing is staged unless the
user asks for it.

    python3 scripts/make_manifest.py --collection 2504 -o config/manifest.tsv
    python3 scripts/make_manifest.py --collection 2504 --samples validation/hprc_overlap.txt \
        -o config/manifest.hprc.tsv

Output is ``sample_id  cram  crai  population``, tab separated, with a header.

Standard library only: this is meant to run on a login node with no environment.
"""

from __future__ import annotations

import argparse
import csv
import sys
import urllib.request
from pathlib import Path

COLLECTIONS = {
    "2504": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
            "1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index",
    "698": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
           "1000G_698_related_high_coverage/1000G_698_related_high_coverage.sequence.index",
}

# The index's own column names, from its '#' header line.
COL_PATH, COL_SAMPLE, COL_POP = "#ENA_FILE_PATH", "SAMPLE_NAME", "POPULATION"


def to_https(ftp_path: str) -> str:
    """ftp://ftp.sra.ebi.ac.uk/vol1/... -> https://ftp.sra.ebi.ac.uk/vol1/...

    Same host, same tree; EBI serves both. HTTPS is what htslib's libcurl
    backend supports and what gets through a firewall that blocks FTP.
    """
    if ftp_path.startswith("ftp://"):
        return "https://" + ftp_path[len("ftp://"):]
    return ftp_path


def read_index(source: str) -> list[dict]:
    """Parse a sequence index from a URL or a local path.

    ``##`` lines are commentary; the single ``#`` line is the header.
    """
    if source.startswith(("http://", "https://")):
        with urllib.request.urlopen(source, timeout=120) as fh:
            text = fh.read().decode()
    else:
        text = Path(source).read_text()

    lines = [ln for ln in text.splitlines() if not ln.startswith("##")]
    if not lines or not lines[0].startswith("#"):
        raise SystemExit(f"{source}: no '#' header line found — is this a sequence index?")
    return list(csv.DictReader(lines, delimiter="\t"))


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--collection", default="2504", choices=sorted(COLLECTIONS),
                   help="2504 unrelated panel (default) or the 698 related samples")
    p.add_argument("--index", help="a sequence index URL or local file, overriding --collection")
    p.add_argument("--samples", type=Path,
                   help="file of sample IDs, one per line; keep only these")
    p.add_argument("-o", "--output", type=Path, required=True)
    args = p.parse_args()

    source = args.index or COLLECTIONS[args.collection]
    rows = read_index(source)

    keep = None
    if args.samples:
        keep = {ln.strip() for ln in args.samples.read_text().splitlines() if ln.strip()}

    seen: dict[str, tuple[str, str]] = {}
    for r in rows:
        path = r.get(COL_PATH, "")
        sample = r.get(COL_SAMPLE, "")
        if not sample or not path.endswith(".cram"):
            continue
        if keep is not None and sample not in keep:
            continue
        # The index has one row per run; high-coverage samples have a single
        # final CRAM, but assert rather than assume, because silently taking the
        # last row would pick an arbitrary one if that ever stopped being true.
        url = to_https(path)
        if sample in seen and seen[sample][0] != url:
            raise SystemExit(
                f"{sample}: two different CRAMs in the index\n"
                f"  {seen[sample][0]}\n  {url}\n"
                "This tool assumes one final CRAM per sample; it is not true here."
            )
        seen[sample] = (url, r.get(COL_POP, ""))

    if keep is not None:
        missing = sorted(keep - set(seen))
        if missing:
            print(f"warning: {len(missing)} requested samples absent from {args.collection}: "
                  f"{', '.join(missing[:8])}{' ...' if len(missing) > 8 else ''}",
                  file=sys.stderr)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample_id", "cram", "crai", "population"])
        for sample in sorted(seen):
            url, pop = seen[sample]
            w.writerow([sample, url, url + ".crai", pop])

    print(f"wrote {len(seen)} samples -> {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
