#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Pick, for every gene, the transcript with the longest CDS (≥ 100 nt on each end,
CDS length ≡ 0 mod 3).  Write those sequences to longest_cDNA.fa, create Bowtie
and HISAT2 indices, and dump a .chrom length file.

"""

from __future__ import annotations
import argparse
import os
import subprocess
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# --------------------------------------------------------------------------- #
#  Core function                                                              
# --------------------------------------------------------------------------- #

def longest_gene(species: str,
                 root: Path = Path("/media/hp/disk4/wzy/reference/hg38"),
                 hisat_dir: Path = Path("/home/hp/miniconda3/bin/hisat2"),
                 bowtie_dir: Path = Path("/home/hp/miniconda3/bin/bowtie"),
                 threads: int = 16) -> None:
    """
    Re‑implementation of longest.gene() from the R script.
    """
    print(f"\n=== processing {species} ===")

    ref_dir = root / species / f"{species}_ref"
    cds_fa  = ref_dir / "CDS_DNA.fa"
    size_tsv = ref_dir / "CDS_intron_size.txt"

    if not cds_fa.exists() or not size_tsv.exists():
        raise FileNotFoundError(f"Missing CDS or size table in {ref_dir}")

    # --------------------------------------------------------------------- #
    # 1. Load CDS sequences and pre‑filter by (len‑200) > 0 & ≡ 0 (mod 3)    
    # --------------------------------------------------------------------- #
    records: list[SeqRecord] = list(SeqIO.parse(cds_fa, "fasta"))
    keep_ids = {
        rec.id
        for rec in records
        if (len(rec) - 200) > 0 and (len(rec) - 200) % 3 == 0
    }

    # --------------------------------------------------------------------- #
    # 2. Read size table and pick the longest CDS per gene                  
    # --------------------------------------------------------------------- #
    df = pd.read_table(size_tsv, dtype=str)          # keep strings -> later cast
    df["id"] = df["gene_id"] + "_" + df["locus"]
    df = df[df["id"].isin(keep_ids)].copy()

    # make sure numeric columns *are* numeric
    df["CDS"] = pd.to_numeric(df["CDS"])
    df["transcript"] = pd.to_numeric(df["transcript"])

    df_longest = (
        df.sort_values(["gene_id", "CDS", "transcript"],
                       ascending=[True, False, False])
          .drop_duplicates("gene_id")
          .sort_values("gene_id")
    )

    # --------------------------------------------------------------------- #
    # 3. Extract matching sequences, remove any duplicates, sort, write out 
    # --------------------------------------------------------------------- #
    rec_map = {rec.id: rec for rec in records}
    sel_records = [
        rec_map[rid] for rid in df_longest["id"] if rid in rec_map
    ]

    out_fa = ref_dir / "longest_cDNA.fa"
    SeqIO.write(sorted(sel_records, key=lambda r: r.id), out_fa, "fasta")
    print(f"  → wrote {len(sel_records)} sequences to {out_fa.name}")

    # --------------------------------------------------------------------- #
    # 4. Build HISAT2 & Bowtie indices                                      
    # --------------------------------------------------------------------- #
    os.environ["PATH"] += f":{hisat_dir}:{bowtie_dir}"

    hisat_prefix  = ref_dir / "longest_hisat"
    bowtie_prefix = ref_dir / "longest"

    subprocess.run(
        ["hisat2-build", "-p", str(threads), out_fa, hisat_prefix],
        check=True
    )
    subprocess.run(
        ["bowtie-build", "--threads", str(threads), out_fa, bowtie_prefix],
        check=True
    )
    print("  → indices built")

    # --------------------------------------------------------------------- #
    # 5. Write “chrom” file (id & CDS+200 length)                           
    # --------------------------------------------------------------------- #
    chrom_df = df_longest[["id", "CDS"]].copy()
    chrom_df["CDS"] = chrom_df["CDS"] + 200
    chrom_path = ref_dir / "longest_cDNA.chrom"
    chrom_df.to_csv(chrom_path, sep="\t", header=False, index=False)
    print(f"  → {chrom_path.name} written")


# --------------------------------------------------------------------------- #
#  CLI                                                                       
# --------------------------------------------------------------------------- #

def main() -> None:
    parser = argparse.ArgumentParser(
        description="creat the longest_cDNA.fa containing the longest isoforms in this file and make the bowtie/hisat2 index for the longest_cDNA"
    )
    parser.add_argument(
        "species", nargs="+",
        help="One or more species folder names (e.g. C_elegans_Ensl_WBcel235)"
    )
    args = parser.parse_args()

    for sp in args.species:
        longest_gene(sp)


if __name__ == "__main__":
    main()
