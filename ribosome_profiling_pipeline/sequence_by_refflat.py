#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Extract CDS‑flanking sequences from a UCSC refFlat annotation, translate them,
count codon usage, and write several helper files.

Direct port of the original Perl script by wz‑y (2025‑07‑07).
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict

FLANK = 100

CODON2AA = {
    # UPPER‑CASE DNA codons → one‑letter AA codes
    "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
    "TTC": "F", "TTT": "F", "TTA": "L", "TTG": "L",
    "TAC": "Y", "TAT": "Y", "TAA": "*", "TAG": "*",
    "TGC": "C", "TGT": "C", "TGA": "*", "TGG": "W",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CAC": "H", "CAT": "H", "CAA": "Q", "CAG": "Q",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
    "ATA": "I", "ATC": "I", "ATT": "I", "ATG": "M",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
    "AAC": "N", "AAT": "N", "AAA": "K", "AAG": "K",
    "AGC": "S", "AGT": "S", "AGA": "R", "AGG": "R",
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
    "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
}

# ---------- helper functions -------------------------------------------------


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of an (upper‑case) DNA sequence."""
    trans = str.maketrans("ATCG", "TAGC")
    return seq.translate(trans)[::-1]


def translate_codon(cod: str) -> str:
    """Translate a single codon; unknown/invalid codons → 'X'."""
    aa = CODON2AA.get(cod.upper())
    if aa is None:
        sys.stderr.write(f"warning: {cod} is not a valid codon\n")
        return "X"
    return aa


def load_genome_fasta(fasta_path: Path) -> dict[str, str]:
    """Load a multi‑FASTA file into a dict {header_without_gt: sequence (UPPER)}."""
    genome: dict[str, list[str]] = defaultdict(list)
    curr = None
    with fasta_path.open() as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                curr = line.split()[0][1:]
            elif curr:
                genome[curr].append(line.upper())
    return {k: "".join(v) for k, v in genome.items()}


# ---------- core logic --------------------------------------------------------


def genome_to_rna(species: str, root: Path) -> None:
    """
    Re‑implement the original Perl sub genomeToRNA().
    Writes the same four output files into <species>/<species>_ref/.
    """
    ref_dir = root / species / f"{species}_ref"
    ref_dir.mkdir(parents=True, exist_ok=True)

    genome_fa = root / "reference" / species / "genome.fa"
    refflat   = root / "reference" / species / "refFlat.txt"

    genome = load_genome_fasta(genome_fa)
    print(f"Finished loading {species} genome ({len(genome):,} contigs)")

    codon_counts: dict[str, int] = defaultdict(int)

    with (
        refflat.open()                         as ref,
        (ref_dir / "CDS_DNA.fa").open("w")     as dna_out,
        (ref_dir / "CDS_pep.fa").open("w")     as pep_out,
        (ref_dir / "translate_warning.txt").open("w") as warn_out,
        (ref_dir / "CDS_intron_size.txt").open("w")   as size_out,
    ):
        size_out.write("gene_id\tlocus\ttranscript\tCDS\tCDS_intron\n")

        for line in ref:
            a = line.rstrip().split("\t")
            chrom = a[2]
            if chrom not in genome:          # genome fragment missing
                continue

            exon_count = int(a[8])
            ins  = [int(x) for x in a[9].rstrip(",").split(",")]
            ine  = [int(x) for x in a[10].rstrip(",").split(",")]
            cds_start, cds_end = int(a[6]), int(a[7])
            tx_start, tx_end   = int(a[4]), int(a[5])
            strand             = a[3]
            gene_id            = a[0]
            transcript_id      = a[1]

            # -----------------------------------------------------------------
            # build full exon coordinate list (0‑based half‑open like UCSC)
            cds_coords = [pos
                          for i in range(exon_count)
                          for pos in range(ins[i], ine[i])]

            if cds_end - cds_start <= 0:
                # no CDS annotated → skip
                continue

            full_len = len(cds_coords)

            # ---------- left flank (3'‑UTR or upstream) ----------------------
            utrL = set(range(tx_start, cds_start))
            left_inter = sorted(utrL & set(cds_coords))
            cds_coords = sorted(set(cds_coords))

            if cds_start == tx_start:
                # No 3′‑UTR annotated → prepend upstream FLANK nt
                if cds_start - FLANK < 0:
                    continue
                cds_coords = list(range(cds_start - FLANK, cds_start)) + cds_coords
            else:
                # Trim or extend so that we keep exactly FLANK nt of UTR
                diff = len(left_inter) - FLANK
                if diff >= 0:
                    cds_coords = cds_coords[diff:]
                else:
                    cds_coords = list(range(tx_start + diff, tx_start)) + cds_coords

            # ---------- right flank (5'‑UTR or downstream) -------------------
            utrR = set(range(cds_end, tx_end))
            right_inter = sorted(utrR & set(cds_coords))

            if cds_end == tx_end:
                # No 5′‑UTR annotated → append downstream FLANK nt
                if cds_end + FLANK > len(genome[chrom]):
                    continue
                cds_coords += list(range(cds_end, cds_end + FLANK))
            else:
                diff = len(right_inter) - FLANK
                if diff >= 0:
                    cds_coords = cds_coords[:-diff]
                else:
                    cds_coords += list(range(tx_end, tx_end + abs(diff)))

            cds_coords = sorted(set(cds_coords))
            cds_core_len = len(cds_coords) - 2 * FLANK
            if cds_core_len <= 0:
                continue

            intron_len = full_len - cds_core_len
            size_out.write(f"{gene_id}\t{strand}\t{full_len}"
                           f"\t{cds_core_len}\t{intron_len}\n")

            # ---------- fetch sequence ---------------------------------------
            seq = "".join(genome[chrom][i] for i in cds_coords)
            if strand == "-":
                seq = reverse_complement(seq)

            # ---------- translate -------------------------------------------
            prot = "".join(
                translate_codon(seq[i:i + 3])
                for i in range(FLANK, len(seq) - FLANK, 3)
            )

            # ---------- classification --------------------------------------
            header = f">{gene_id}_{transcript_id}_{strand}"
            if "X" not in prot and prot.count("*") <= 1:   # no internal stop/X
                dna_out.write(f"{header}\n{seq}\n")
                pep_out.write(f"{header}\n{prot}\n")
                for i in range(FLANK, len(seq) - FLANK, 3):
                    cod = seq[i:i + 3]
                    codon_counts[cod] += 1
            else:
                warn_out.write(f"{header}\n{prot}\n{seq}\n")

    # ---------- write codon frequency table ----------------------------------
    freq_file = ref_dir / "codon_frequency.txt"
    with freq_file.open("w") as fh:
        for cod, cnt in sorted(codon_counts.items()):
            fh.write(f"{cod}\t{translate_codon(cod)}\t{cnt}\n")


# ---------- entry point ------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract sequences from a UCSC refFlat file (Python rewrite)."
    )
    parser.add_argument(
        "-g", "--genome", required=True,
        help="species / folder name (e.g. C_elegans_Ensl_WBcel235)"
    )
    parser.add_argument(
        "-r", "--root", default="/media/hp/disk4/wzy",
        help="root directory that contains 'reference/' (default: %(default)s)"
    )
    args = parser.parse_args()

    root = Path(args.root).resolve()
    species = args.genome

    genome_to_rna(species, root)


if __name__ == "__main__":
    main()
