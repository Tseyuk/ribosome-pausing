STEP 1 Reference Building
1 processing genome fasta
remove all unwanted charaters in the chromosome header of fasta file.The genome fasta if downloaded from NCBI(https://www.ncbi.nlm.nih.gov/genome/) or UCSC(https://hgdownload.soe.ucsc.edu/downloads.html)

sed -i -e 's/ .*//' genome.fa
2 create refFlat.txt
using UCSC tools gtfToGenePred (http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64)

~/Genomes/tools/gtfToGenePred -genePredExt genes.gtf ref.txt
move the 12th column to first

awk 'BEGIN {FS=OFS="\t"} {$1=$12FS$1;NF--} 1' ref.txt >ref1.txt
remove the last 4 columns

awk 'NF{NF-=4}1' FS="\t" OFS="\t" ref1.txt > refFlat.txt
3 build reference
creat the CDS.fa like this 100nt 5'-UTR + CDS + 100nt 3'-UTR

perl cDNA/sequence_by_refFlat.pl  -g ${name}
creat the longest_cDNA.fa containing the longest isoforms in this file and make the bowtie/hisat2 index for the longest_cDNA

Rscript cDNA/longest_cDNA_index.R ${name}
4 calculate the size of genome
samtools faidx genome.fa
cut -f1,2 genome.fa.fai > sizes.genome
5 calculate CBI
Rscript cDNA/CBI.R ${name}				
perl cDNA/UTR_length.pl -g ${name}
perl cDNA/cal_CBI_CAIavg.pl -sp ${name}	
6 creat index for genome
python ~/Genomes/extract_exons.py genes.gtf >genome.exon
python ~/Genomes/extract_splice_sites.py genes.gtf >genome.ss
hisat2-bulid -p 16 genome.fa --ss ~/Genomes/${name}/Genes/genome.ss --exon ~/Genomes/${name}/Genes/genome.exon genome
7 creat index for rRNA
hisat2-bulid -p 16 ribosomes.fa rDNA
