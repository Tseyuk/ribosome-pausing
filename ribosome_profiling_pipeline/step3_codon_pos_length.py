#statistic of 3-nt periodicity based the length of reads.
import os
from collections import defaultdict
from decimal import Decimal, ROUND_HALF_UP

def round_half_up(n):
    return int(Decimal(n).quantize(0, ROUND_HALF_UP))

# Parameters
up = 40
low = 25

# Create output directory if it doesn't exist
os.makedirs("codon_stat", exist_ok=True)

# Read filenames from 'deadapter.txt'
with open("deadapter.txt", "r") as file:
    for line in file:
        sample = line.strip()
        print(sample)
        process_sample(sample)

def process_sample(sam):
    count = 0
    count_dict = defaultdict(int)
    p0 = defaultdict(int)
    pplus = defaultdict(int)
    pminus = defaultdict(int)
    size_dict = defaultdict(int)
    strand_dict = defaultdict(int)
    tot = {'0': 0, 'pplus': 0, 'pminus': 0}

    try:
        with open(f"{sam}.sort.bed", "r") as hand1:
            for line in hand1:
                line = line.rstrip()
                a1 = line.split('\t')
                
                strand = a1[5]
                strand_dict[strand] += 1

                if strand == '-':
                    continue

                start = int(a1[1])
                end = int(a1[2])
                size = end - start
                size_dict[size] += 1

                if size < low or size > up:
                    continue

                count_dict[size] += 1
                count += 1

                f12 = (start + 15 - 100) % 3

                if f12 == 0:
                    tot['0'] += 1
                    p0[size] += 1
                elif f12 == 1:
                    tot['pplus'] += 1
                    pplus[size] += 1
                else:
                    tot['pminus'] += 1
                    pminus[size] += 1

    except FileNotFoundError:
        print(f"File {sam}.sort.bed not found!")
        return

    # Write per-sample output
    with open(f"codon_stat/{sam}.txt", "w") as hand2:
        for key, val in strand_dict.items():
            hand2.write(f"{key}\t{val}\n")

        for key in sorted(size_dict.keys()):
            hand2.write(f"{key}\t{size_dict[key]}\n")

        hand2.write(f"The mapped reads of {low}-{up} nt are {count}\n")
        hand2.write("length\tnumber_of_reads\t0\tplus\tminus\n")

        for si in range(low, up + 1):
            if count_dict[si] > 0:
                ra12 = round_half_up(p0[si] / count_dict[si] * 100) if count_dict[si] else 0
                ra13 = round_half_up(pplus[si] / count_dict[si] * 100) if count_dict[si] else 0
                ra14 = round_half_up(pminus[si] / count_dict[si] * 100) if count_dict[si] else 0
                if size_dict[si]:
                    hand2.write(f"{si}\t{size_dict[si]}\t{ra12}\t{ra13}\t{ra14}\n")

    # Append total summary
    total_file = "codon_stat/total.txt"
    write_header = not os.path.exists(total_file)

    with open(total_file, "a") as hand3:
        if write_header:
            hand3.write("sample\t0\tplus\tminus\n")
        hand3.write(f"{sam}\t{tot['0']}\t{tot['pplus']}\t{tot['pminus']}\n")

    print(f"{sam} is done")
