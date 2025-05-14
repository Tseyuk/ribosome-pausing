#step4_trimbed for generating _A.bed files
import os
import argparse
from collections import defaultdict

# Argument parsing
parser = argparse.ArgumentParser(description="Trim RPF fragments for codon position analysis.")
parser.add_argument('-t', type=int, default=60, help='Threshold for codon assignment')
parser.add_argument('-sp', type=str, required=True, help='Species/sample name')
parser.add_argument('-up', type=int, default=40, help='Upper bound of fragment size')
parser.add_argument('-down', type=int, default=25, help='Lower bound of fragment size')

args = parser.parse_args()
threshold = args.t
species = args.sp
up = args.up
down = args.down

if not species:
    raise SystemExit("Missing species name. Use -sp <species name>")

os.makedirs('cleavage_preference', exist_ok=True)

# Load transcriptome reference
gene = {}
si = {}

ref_path = f"/disk/wzy/ribo_seq/human_ribo/{species}/longest_cDNA.fa"
try:
    with open(ref_path, "r") as ref:
        ge = None
        for line in ref:
            line = line.strip()
            if line.startswith(">"):
                ge = line[1:]
                gene[ge] = ""
            else:
                gene[ge] += line
    for k in gene:
        si[k] = len(gene[k])
except FileNotFoundError:
    raise SystemExit(f"Reference file not found: {ref_path}")

# Read sample names
with open("deadapter.txt", "r") as f:
    samples = [line.strip() for line in f if line.strip()]

for sample in samples:
    print(f"Processing {sample}")
    def trimbed(sam):
        trim = {}
        count_dict = defaultdict(int)
        pref5 = defaultdict(lambda: defaultdict(int))

        # Step 1: Read codon_stat file
        try:
            with open(f"codon_stat/{sam}.txt", "r") as h0:
                for line in h0:
                    if not line.startswith(tuple(str(x) for x in range(10))):  # skip non-numeric start
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) < 5:
                        continue
                    length = int(parts[0])
                    val0, val1, val2 = int(parts[2]), int(parts[3]), int(parts[4])
                    if val0 >= threshold:
                        trim[length] = 15
                    elif val1 >= threshold:
                        trim[length] = 14
                    elif val2 >= threshold:
                        trim[length] = 16
        except FileNotFoundError:
            print(f"File codon_stat/{sam}.txt not found.")
            return

        # Step 2: Process BED file
        try:
            with open(f"{sam}.sort.bed", "r") as h1, \
                 open(f"{sam}_A.bed", "w") as h3, \
                 open(f"codon_stat/{sam}_trim.txt", "w") as h4, \
                 open(f"cleavage_preference/{sam}.5end.stat.txt", "w") as h5:

                h4.write(f"trimming status >= {threshold}\nlength\tposition\n")
                h5.write("Junction_4nt\tInner\tinframe\tplus1\tminus1\n")

                for line in h1:
                    parts = line.strip().split('\t')
                    if len(parts) < 6:
                        continue
                    chrom, start, end, name, score, strand = parts
                    if strand == '-':
                        continue

                    start, end = int(start), int(end)
                    frag_size = end - start
                    if frag_size < down or frag_size > up:
                        continue
                    if frag_size not in trim:
                        continue

                    count_dict[frag_size] += 1
                    pos = start + trim[frag_size] - 3
                    si_val = pos + 1
                    posA = pos + 3
                    siA = posA + 1

                    h3.write(f"{chrom}\t{posA}\t{siA}\t{name}\t{score}\t{strand}\n")

                    if pos < 100 or (pos + 100 > si.get(chrom, 0)):
                        continue

                    try:
                        junc5 = gene[chrom][start - 2: start + 2]
                    except KeyError:
                        continue

                    frame = "in"
                    if (pos - 100) % 3 == 1:
                        frame = "s1"
                    elif (pos - 100) % 3 == 2:
                        frame = "s2"
                    pref5[junc5][frame] += 1

                for k in sorted(trim.keys()):
                    if k in count_dict:
                        h4.write(f"{k}\t{trim[k]}\t{count_dict[k]}\n")

                for junc4 in pref5:
                    inner = junc4[1:3]
                    h5.write(f"{junc4}\t{inner}")
                    for state in ["in", "s1", "s2"]:
                        h5.write(f"\t{pref5[junc4].get(state, 0)}")
                    h5.write("\n")

        except FileNotFoundError:
            print(f"Missing input file: {sam}.sort.bed")

    trimbed(sample)
