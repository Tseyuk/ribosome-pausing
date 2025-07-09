# Ribosome-pausing
This script is designed to quickly generate a file recording all the signals near to the central 3 codons flanked by 10 codons by using sliding window. It can also be running at  underpowered or RAM-limited machine. 

## Ribosome-pausing/

```shell
├── Script in python for analyzing ribosome pausing/
│   ├── test_A.bedgraph           # bedgraph  demo file for testing
│   ├── test_longest_cDNA.txt     # longest_cDNA demo file as reference
│   ├── tri_PPP_pausing.png       # image of a typical pausing position
│   ├── v3_finnal_ave_0frame_normed_dicodon_pausing.py   # Main script to run the pipeline
│   └── Readme.txt               
│
├── ribosome_profiling_pipeline/  # to generate .bedgraph file
│
├── README.md                     # Documentation (this file)
```

## Usage

python subsampling.py v3_finnal_ave_0frame_normed_dicodon_pausing.py -i  -r 

-i, input .bedgraph file

-r, input longest_cDNA.txt file

## Example:

python v3_finnal_ave_0frame_normed_dicodon_pausing.py -i .\test_A.bedgraph -r .\ .\test_longest_cDNA.txt

NOTE:test_A.bedgraph is a file that records ribosome footprints
test_longest_cDNA.txt is a file that contains the lonest cDNA for a gene
tri_PPP_pausing.png is a figure shows that there is a strong pausing when translating three Prolines in a cell line that the eIF5A is knockdown, and eIF5A is crucial for PPP decoding, so there is a strong stalling compared with the WT.
