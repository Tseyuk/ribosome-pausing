## Installation

When all steps in ribosome_profiling_pipeline\ have been conducted, a bedgraph file will be created, then 

At terminal, use

python v3_finnal_ave_0frame_normed_dicodon_pausing.py -i .\test_A.bedgraph -r .\ .\test_longest_cDNA.txt

NOTE:test_A.bedgraph is a file that records ribosome footprints
test_longest_cDNA.txt is a file that contains the lonest cDNA for a gene
tri_PPP_pausing.png is a figure shows that there is a strong pausing when translating three Prolines in a cell line that the eIF5A is knockdown, and eIF5A is crucial for PPP decoding, so there is a strong stalling compared with the WT.
