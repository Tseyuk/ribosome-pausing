#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
bamtobed.py – Python 3 rewrite of the original BAM‑to‑BED converter.

* For every *sample* listed in  <root>/b_align/deadapter.txt
  ‑ Read `<root>/b_align/codon_stat/<sample>.txt` to decide the offset
    (+15 / +14 / +16) for each read length whose in‑frame counts exceed
    the `--periodicity` threshold.
  ‑ Stream‑parse `<root>/align2/<sample>.bam` (via `samtools view`),
    convert each read’s P‑site to a 1‑bp BED record at the appropriate
    offset, and write `<root>/align2/<sample>.bed`.
  ‑ Mis‑flagged reads are copied to
    `<root>/b_align/codon_stat/<sample>_other.txt`.
"""

from __future__ import annotations
import argparse
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List

# --------------------------------------------------------------------------- #
#  Helpers                                                                     #
# --------------------------------------------------------------------------- #

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")   # CIGAR element parser


def mapping_length(cigar: str) -> int:
    """Sum M and D lengths in a CIGAR string (same as Perl regexp)."""
    return sum(int(n) for n, op in CIGAR_RE.findall(cigar) if op in "MD")


def iter_bam(path: Path) -> List[str]:
    """Yield SAM lines from `samtools view <bam>`."""
    proc = subprocess.Popen(
        ["samtools", "view", str(path)],
        stdout=subprocess.PIPE, text=True, bufsize=1
    )
    assert proc.stdout is not None           # for mypy
    for line in proc.stdout:
        yield line.rstrip("\n")
    proc.stdout.close()
    if proc.wait():
        sys.exit(f"samtools view failed on {path}")


# --------------------------------------------------------------------------- #
#  Per‑sample workflow                                                         #
# --------------------------------------------------------------------------- #

def offset_table(stat_file: Path, threshold: int) -> Dict[int, int]:
    """
    Read <sample>.txt from codon_stat and return {read_len: offset}.
    Offset is 15 if in‑frame counts ≥ threshold,
    14 if +1 frame ≥ threshold,
    16 if –1 frame ≥ threshold.
    """
    tbl: Dict[int, int] = {}
    with stat_file.open() as fh:
        for ln in fh:
            ln = ln.rstrip()
            if not re.match(r"^\d{2}\t\d+\t", ln):
                continue
            length_s, num, in0_s, plus_s, minus_s = ln.split("\t")
            length = int(length_s)
            in0   = int(in0_s)
            plus  = int(plus_s)
            minus = int(minus_s)
            if in0 >= threshold:
                tbl[length] = 15
            elif plus >= threshold:
                tbl[length] = 14
            elif minus >= threshold:
                tbl[length] = 16
    return tbl


def process_sample(sample: str, root: Path, periodicity: int) -> None:
    print(f"\nprocessing {sample} …")

    stat_file = root / "b_align" / "codon_stat" / f"{sample}.txt"
    bam_file  = root / "align2"   / f"{sample}.bam"
    bed_file  = root / "align2"   / f"{sample}.bed"
    other_file = stat_file.with_name(f"{sample}_other.txt")

    if not stat_file.exists() or not bam_file.exists():
        sys.exit(f"Missing codon_stat or BAM for sample {sample}")

    trim = offset_table(stat_file, periodicity)
    if not trim:
        print(f"  (no read lengths pass threshold {periodicity})")
        return
    print(f"  read‑length → offset table: {trim}")

    # Open outputs
    bed_out   = bed_file.open("w")
    other_out = other_file.open("a")

    for sam in iter_bam(bam_file):
        fields = sam.split("\t")
        qname, flag_s, rname, pos_s, mapq, cigar = fields[:6]
        flag = int(flag_s)
        if flag & 4:        # unmapped
            continue

        strand = "+" if flag & 16 == 0 else "-"
        bam_start = int(pos_s) - 1
        read_len = mapping_length(cigar)

        off = trim.get(read_len)
        if off is None:
            continue

        # Parse full CIGAR for splice-aware offset
        elements = [(int(n), op) for n, op in CIGAR_RE.findall(cigar)]
        splice = 0
        seen   = 0

        if strand == "+":
            for n, op in elements:
                seen += n
                if op in "MD":
                    splice += n
                if splice > 15:
                    skip = off - (splice - n)
                    bed_start = bam_start + (seen - n) + skip
                    bed_out.write(
                        f"{rname}\t{bed_start}\t{bed_start+1}\t{qname}\t{mapq}\t+\n"
                    )
                    break
        else:   # strand == "-"
            bam_end = bam_start + sum(n for n, _ in elements)
            for n, op in reversed(elements):
                seen += n
                if op in "MD":
                    splice += n
                if splice > 15:
                    skip = off - (splice - n)
                    bed_end   = bam_end - (seen - n) - skip
                    bed_start = bed_end - 1
                    bed_out.write(
                        f"{rname}\t{bed_start}\t{bed_end}\t{qname}\t{mapq}\t-\n"
                    )
                    break
        # unexpected flags
        if flag & 3840 not in (0, 16, 256, 272):   # primary/secondary ±strand
            other_out.write(sam + "\n")

    bed_out.close()
    other_out.close()
    print(f"  BED written ➜ {bed_file}")


# --------------------------------------------------------------------------- #
#  CLI                                                                         #
# --------------------------------------------------------------------------- #

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract P‑site BED positions from BAMs (ribosome profiling)"
    )
    parser.add_argument("-d", "--doc", required=True,
                        help="root directory containing align2/ and b_align/")
    parser.add_argument("-p", "--periodicity", type=int, required=True,
                        help="in‑frame read count threshold")
    args = parser.parse_args()

    root = Path(args.doc).resolve()
    deadapter = root / "b_align" / "deadapter.txt"
    if not deadapter.exists():
        sys.exit(f"{deadapter} not found")

    with deadapter.open() as fh:
        samples = [ln.strip() for ln in fh if ln.strip()]

    for s in samples:
        process_sample(s, root, args.periodicity)
        print(f"=================== finished {s} ===================")


if __name__ == "__main__":
    main()

