#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
step1_trim3.py  –  Python 3 rewrite of step1_trim3.pl (2020‑09‑14).

* Reads a `filelist.txt` (one sample id per line, no extension).
* For every `<sample>.fq.gz`:
    – Trim the first occurrence of the adapter (or its 1‑5 nt‑shorter head)
      from the *3′ end* of each read.
    – Keep reads longer than `--min` (default 20 nt).
    – Crop quality strings to the new length.
    – Write:
        •  <sample>/<sample>_<min>.txt        (trimmed sequences)
        •  <sample>/<sample>_<min>_uni.txt    (unique seqs, FASTA)
        •  <sample>_<min>.fastq               (trimmed FASTQ)
        •  <sample>/sum_<sample>.txt          (stats)
"""

from __future__ import annotations
import argparse
import gzip
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Tuple

# --------------------------------------------------------------------------- #
#  Helper functions                                                           
# --------------------------------------------------------------------------- #

def trim_adapter(seq: str, adapter: str) -> str:
    """
    Return `seq` truncated at the first occurrence of `adapter` (or adapter
    missing the last 1‑5 nt).  If none found, return the original sequence.
    """
    pos = seq.find(adapter)
    if pos != -1:
        return seq[:pos]

    # allow 1‑5 nt 3′ mismatches
    for n_miss in range(1, 6):
        prefix = adapter[:-n_miss]
        if len(prefix) < 3:       # nothing sensible left to match
            break
        pos = seq.find(prefix)
        if pos != -1:
            return seq[:pos]
    return seq


def iter_fastq_gz(path: Path) -> Iterable[Tuple[str, str, str, str]]:
    """Yield 4‑tuples (header, seq, sep, qual) from a gzipped FASTQ file."""
    with gzip.open(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline()
            sep = fh.readline()
            q = fh.readline()
            if not q:
                raise ValueError(f"Broken FASTQ: {path}")
            yield h.rstrip(), s.rstrip(), sep.rstrip(), q.rstrip()


# --------------------------------------------------------------------------- #
#  Per‑sample processing                                                      
# --------------------------------------------------------------------------- #

def process_sample(sample: str, adapter: str, min_len: int, id_start: int = 10_000_000) -> None:
    print(f"\nprocessing {sample}: trimming 3′ adapter")

    fq_in  = Path(f"{sample}.fq.gz")
    if not fq_in.exists():
        sys.exit(f"Input file {fq_in} not found")

    # — outputs & bookkeeping — ------------------------------------------------
    out_dir = Path(sample)
    out_dir.mkdir(exist_ok=True)

    seq_txt   = out_dir / f"{sample}_{min_len}.txt"
    fq_out    = Path(f"{sample}_{min_len}.fastq")
    uni_fa    = out_dir / f"{sample}_{min_len}_uni.txt"
    summary   = out_dir / f"sum_{sample}.txt"

    total_kept  = 0
    first_base  = Counter()
    size_hist   = Counter()
    seq_counts  = Counter()

    next_id = id_start

    with seq_txt.open("w") as seq_handle, fq_out.open("w") as fq_handle:
        for head, seq, sep, qual in iter_fastq_gz(fq_in):
            trimmed = trim_adapter(seq, adapter)
            if len(trimmed) <= min_len:
                continue

            total_kept += 1
            first_base.update(trimmed[0])
            size_hist.update([len(trimmed)])
            seq_counts.update([trimmed])

            seq_handle.write(f"{trimmed}\n")

            next_id += 1
            fq_handle.write(f"{head}\n{trimmed}\n+{next_id}\n{qual[:len(trimmed)]}\n")

    # — write unique fasta — ---------------------------------------------------
    with uni_fa.open("w") as fh:
        unique_id = id_start
        for seq, cnt in seq_counts.most_common():
            unique_id += 1
            fh.write(f">{unique_id}_x{cnt}\n{seq}\n")

    # — write summary — --------------------------------------------------------
    with summary.open("w") as fh:
        fh.write(f"{sample}\n\n")
        fh.write(f"Total number of {min_len} nt or longer reads\t{total_kept}\n")
        fh.write(f"Total unique species of {min_len} nt or longer\t{len(seq_counts)}\n")

        fh.write("\nstart:\n")
        for base, cnt in first_base.items():
            fh.write(f"{base}\t{cnt}\n")

        fh.write("\nsize:\n")
        for size in sorted(size_hist):
            fh.write(f"{size}\t{size_hist[size]}\n")

    print(f"  kept {total_kept:,} reads  ➜  {fq_out}")
    print(f"  unique seqs: {len(seq_counts):,}  ➜  {uni_fa}")
    print(f"  summary      ➜  {summary}")


# --------------------------------------------------------------------------- #
#  Main                                                                       
# --------------------------------------------------------------------------- #

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Trim 3′ adapters from .fq.gz files listed in filelist.txt"
    )
    parser.add_argument("-ad", "--adapter", required=True,
                        help="adapter sequence to trim (exact, 3′ linker)")
    parser.add_argument("-m", "--min", type=int, default=20,
                        help="minimum read length to keep (default: 20)")
    args = parser.parse_args()

    # read sample list
    if not Path("filelist.txt").exists():
        sys.exit("filelist.txt not found")
    with open("filelist.txt") as fl:
        samples = [ln.strip() for ln in fl if ln.strip()]

    if not samples:
        sys.exit("filelist.txt is empty")

    print("This script is updated on 20200914 by suang to remove the 3' linker\n")
    for s in samples:
        process_sample(s, args.adapter, args.min)


if __name__ == "__main__":
    main()
