"""IUPAC-aware CDS extraction and translation.

Vendored from `lilr-genotyper`'s ``lilr_genotyper/sequence_builder.py`` and kept
close to the original: the translation and CDS-slicing logic is correct and
well-exercised, and diverging from it would make allele sequences from the two
pipelines incomparable — which is the whole point of validating one against the
other.

Only the parts this pipeline reaches are live: :func:`extract_sequences`,
:func:`_translate` and :func:`_revcomp`. ``MIN_CONS_DEPTH`` and
:func:`build_consensus` belong to the upstream pileup path and are *not* used
here; masking is driven by :mod:`lilrwgs.depth_model` instead, which is the
substantive difference between the two pipelines and must not be quietly
reintroduced as a constant.
"""

ORIGINAL_DOCSTRING = """
Build gDNA / gDNA_CDS / cDNA / protein consensus sequences from pileup data,
and assign provisional LILR nomenclature.

Nomenclature format (analogous to HLA/KIR):
  GENE*XXX YY ZZ
  XXX : protein-level allele number   (001–999)
  YY  : synonymous CDS differences    (01–99)
  ZZ  : non-coding / intronic / flank (01–99)

All fields are zero-padded; the full allele designator is GENE*XXXYYZZ.
"""

import json
import logging
import re
import subprocess
import tempfile
import os
from pathlib import Path
from dataclasses import dataclass, field

log = logging.getLogger(__name__)

MIN_CONS_DEPTH  = 5    # min depth to call a consensus base; else mask as N
LOW_CONF_DEPTH  = 10   # depth threshold below which call is flagged low-confidence
NOVEL_MIN_DEPTH = 10   # minimum total depth at a position to call novel
HET_THRESHOLD   = 0.25 # minor-allele fraction above this → heterozygous

IUPAC = {
    frozenset('AC'): 'M', frozenset('AG'): 'R', frozenset('AT'): 'W',
    frozenset('CG'): 'S', frozenset('CT'): 'Y', frozenset('GT'): 'K',
    frozenset('ACG'): 'V', frozenset('ACT'): 'H', frozenset('AGT'): 'D',
    frozenset('CGT'): 'B', frozenset('ACGT'): 'N',
}

# Reverse mapping: IUPAC ambiguity code → set of possible bases (for translation)
IUPAC_BASES = {
    'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'},
    'M': {'A','C'}, 'R': {'A','G'}, 'W': {'A','T'},
    'S': {'C','G'}, 'Y': {'C','T'}, 'K': {'G','T'},
    'V': {'A','C','G'}, 'H': {'A','C','T'},
    'D': {'A','G','T'}, 'B': {'C','G','T'},
    'N': {'A','C','G','T'},
}

CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}


@dataclass
class SequenceResult:
    locus: str
    sample_id: str
    gdna: str = ''          # full reference-length sequence (flanks + gene)
    gdna_cds: str = ''      # from mRNA start to mRNA end (includes introns)
    cdna: str = ''          # spliced exons only
    protein: str = ''       # translated cDNA
    novel_snps: list = field(default_factory=list)  # list of (pos, ref_base, obs_base, depth)
    allele_name: str = ''   # provisional name e.g. LILRA1*00101 01


def load_exon_coords(resources_dir: Path) -> dict:
    path = resources_dir / 'exon_coords.json'
    with open(path) as f:
        return json.load(f)


def _get_ref_sequence(locus: str, resources_dir: Path) -> str:
    """Read the locus reference sequence (single contig)."""
    fa = resources_dir / 'references' / f'{locus}_named.fa'
    seq = ''
    with open(fa) as f:
        for line in f:
            if not line.startswith('>'):
                seq += line.strip()
    return seq


def build_consensus(ref_seq: str, pileup: dict,
                    cn: int = 2,
                    min_depth: int = MIN_CONS_DEPTH,
                    low_conf_depth: int = LOW_CONF_DEPTH,
                    trusted_end: int = None):
    """
    Build a CN-aware consensus sequence over the full reference.

    Positions with depth < min_depth (default 5) are masked as N.
    Positions with min_depth ≤ depth < low_conf_depth (default 10) are called
    but counted as low-confidence (returned via the n_low_confidence value).
    Positions with depth ≥ low_conf_depth are full-confidence calls.

    Returns (consensus_str, n_low_confidence_positions).
    """
    seq = list(ref_seq)
    n_low_conf = 0
    for pos, counts in pileup.items():
        if trusted_end is not None and pos >= trusted_end:
            continue
        if not counts:
            seq[pos] = 'N'
            continue
        total = sum(counts.values())
        if total < min_depth:
            seq[pos] = 'N'
            continue
        if total < low_conf_depth:
            n_low_conf += 1

        if cn == 1:
            seq[pos] = max(counts, key=counts.get)
            continue

        het_bases = frozenset(b for b, c in counts.items()
                              if c / total >= HET_THRESHOLD)
        if len(het_bases) == 0:
            seq[pos] = max(counts, key=counts.get)
        elif len(het_bases) == 1:
            seq[pos] = next(iter(het_bases))
        else:
            seq[pos] = IUPAC.get(het_bases, 'N')
    return ''.join(seq), n_low_conf


_RC = str.maketrans(
    'ACGTacgtMRWSYKVHDBNmrwsykvhdbn',
    'TGCAtgcaKYWSRMBDHVNkywsrmbdhvn',
)
# IUPAC complements: M↔K, R↔Y, W↔W, S↔S, V↔B, H↔D, N↔N

def _revcomp(seq: str) -> str:
    return seq.translate(_RC)[::-1]


def extract_sequences(consensus: str, coords: dict) -> tuple:
    """
    Return (gdna_cds, cdna, protein) from consensus sequence.
    coords: dict with keys mrna_start, mrna_end, strand, exons (sorted by coord).
    """
    mrna_s = coords['mrna_start']
    mrna_e = coords['mrna_end']
    strand = coords.get('strand', '+')
    exons  = coords['exons']   # sorted low→high

    gdna_cds = consensus[mrna_s:mrna_e]
    if strand == '-':
        gdna_cds = _revcomp(gdna_cds)

    if strand == '+':
        cdna = ''.join(consensus[s:e] for s, e in exons)
    else:
        # Minus strand: exons in reverse coordinate order, each rev-complemented
        cdna = ''.join(_revcomp(consensus[s:e]) for s, e in reversed(exons))

    protein = _translate(cdna)
    return gdna_cds, cdna, protein


def _translate(cdna: str) -> str:
    """
    Translate cDNA to protein, handling IUPAC ambiguity and N.

    For each codon, enumerate all concrete codons implied by the IUPAC bases.
    If they all encode the same amino acid (synonymous under ambiguity) the
    codon is translated unambiguously; if they differ, the codon translates
    to X (unknown).  Stops at first stop codon.
    """
    aa = []
    for i in range(0, len(cdna) - 2, 3):
        codon = cdna[i:i+3].upper()
        # Expand IUPAC: get set of possible bases at each codon position
        opts = [IUPAC_BASES.get(b, set()) for b in codon]
        if any(not o for o in opts):
            aa.append('X')
            continue
        # Enumerate all concrete codons
        possible_aa = set()
        for b1 in opts[0]:
            for b2 in opts[1]:
                for b3 in opts[2]:
                    possible_aa.add(CODON_TABLE.get(b1+b2+b3, 'X'))
        if len(possible_aa) == 1:
            aa_char = next(iter(possible_aa))
            if aa_char == '*':
                break
            aa.append(aa_char)
        else:
            # Ambiguous codon — output X unless all are stop (rare)
            if possible_aa == {'*'}:
                break
            aa.append('X')
    return ''.join(aa)


def find_novel_snps(pileup: dict, ref_seq: str,
                    allele1_bases: dict, allele2_bases: dict,
                    min_depth: int = NOVEL_MIN_DEPTH,
                    trusted_end: int = None) -> list:
    """
    Identify positions with depth ≥ min_depth where the consensus (major) base
    differs from both called catalog alleles at that position.

    allele1_bases / allele2_bases: dict[0-based pos] = expected base char
    Returns list of (pos, ref_base, novel_base, read_count).
    """
    novel = []
    for pos, counts in pileup.items():
        if trusted_end is not None and pos >= trusted_end:
            continue
        if not counts:
            continue
        depth = sum(counts.values())
        if depth < min_depth:
            continue
        ref_b  = ref_seq[pos].upper() if pos < len(ref_seq) else 'N'
        exp_b1 = allele1_bases.get(pos, ref_b).upper()
        exp_b2 = allele2_bases.get(pos, ref_b).upper()
        expected = {exp_b1, exp_b2}

        major = max(counts, key=counts.get)
        if major.upper() not in expected:
            novel.append((pos, ref_b, major.upper(), counts[major]))
    return novel


# ---------------------------------------------------------------------------
# Provisional nomenclature
# ---------------------------------------------------------------------------

def _equivalent(a: str, b: str, unknown_chars: frozenset) -> bool:
    """
    Two sequences represent the same allele if:
      - Same length
      - At every position, either side is unknown (N for DNA, X for protein),
        OR the bases agree.
    """
    if len(a) != len(b):
        return False
    for ca, cb in zip(a, b):
        if ca in unknown_chars or cb in unknown_chars:
            continue
        if ca != cb:
            return False
    return True


_DNA_UNKNOWN = frozenset('N')
_PROT_UNKNOWN = frozenset('X')


class NomenclatureRegistry:
    """
    Assigns provisional allele names using a 7-digit code GENE*XXXYYZZ.
      XXX = protein-level index
      YY  = synonymous cDNA variant within the same protein group
      ZZ  = non-coding/flank variant within the same protein+cDNA group

    Uses equivalence-based matching so that two sequences with different
    N/X masking patterns but no conflicting positions are treated as the
    same allele. Tracks per-allele sample counts for recurrence filtering.
    """

    def __init__(self, registry_path: Path):
        self.path = registry_path
        if registry_path.exists():
            with open(registry_path) as f:
                self._data = json.load(f)
        else:
            self._data = {}

    def _locus_data(self, locus: str) -> dict:
        if locus not in self._data:
            self._data[locus] = {
                'proteins': [],  # list of [idx, sequence]
                'cdnas':    [],  # list of [prot_idx, cdna_idx, sequence]
                'gdnas':    [],  # list of [prot_idx, cdna_idx, gdna_idx, sequence]
                'counts':   {},  # allele name → sample count
            }
        return self._data[locus]

    def _find_or_add_protein(self, ld: dict, protein: str) -> int:
        for idx, seq in ld['proteins']:
            if _equivalent(seq, protein, _PROT_UNKNOWN):
                return idx
        new_idx = len(ld['proteins']) + 1
        ld['proteins'].append([new_idx, protein])
        return new_idx

    def _find_or_add_cdna(self, ld: dict, prot_idx: int, cdna: str) -> int:
        existing = [e for e in ld['cdnas'] if e[0] == prot_idx]
        for _, idx, seq in existing:
            if _equivalent(seq, cdna, _DNA_UNKNOWN):
                return idx
        new_idx = len(existing) + 1
        ld['cdnas'].append([prot_idx, new_idx, cdna])
        return new_idx

    def _find_or_add_gdna(self, ld: dict, prot_idx: int, cdna_idx: int,
                          gdna: str) -> int:
        existing = [e for e in ld['gdnas']
                    if e[0] == prot_idx and e[1] == cdna_idx]
        for _, _, idx, seq in existing:
            if _equivalent(seq, gdna, _DNA_UNKNOWN):
                return idx
        new_idx = len(existing) + 1
        ld['gdnas'].append([prot_idx, cdna_idx, new_idx, gdna])
        return new_idx

    def assign(self, locus: str, protein: str, cdna: str, gdna_full: str) -> str:
        """
        Return the allele name for this sequence combination, creating a new
        entry if no equivalent is already registered.  Increments sample count.
        """
        ld = self._locus_data(locus)
        prot_idx = self._find_or_add_protein(ld, protein)
        cdna_idx = self._find_or_add_cdna(ld, prot_idx, cdna)
        gdna_idx = self._find_or_add_gdna(ld, prot_idx, cdna_idx, gdna_full)
        name = f'{locus}*{prot_idx:03d}{cdna_idx:02d}{gdna_idx:02d}'
        ld['counts'][name] = ld['counts'].get(name, 0) + 1
        return name

    def sample_count(self, locus: str, name: str) -> int:
        return self._locus_data(locus)['counts'].get(name, 0)

    def save(self):
        with open(self.path, 'w') as f:
            json.dump(self._data, f)
