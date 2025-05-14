Ribo-seq_process
species="humanref"
doc="/disk/ribo_seq/xxx"
rna="yes"
position="A"
Genome="/disk/ribo_seq/humanref"
sort="longest"   #CDS or longest
sort_file="longest_cDNA"   
cd ${doc}
mkdir b_align
mkdir de_rRNA
cat deadapter.txt | while read line
do
        if [ ${rna} = "yes" ];then
                #remove rRNA
                echo "=============== processing $line discard rRNA ================="
                echo "=============== processing $line =================">> de_rRNA/rRNA_stats.txt
                bowtie -p 56 -v 3 -5 3 --un=${line}_derRNA.fq ${Genome}/ribosome/rDNA ${line}_trimmed.fq.gz 2>>  de_rRNA/rRNA_stats.txt > de_rRNA/${line}_rRNAalign.aln  #note ,sometimes we need "-3 4",and another we may don't need . responding to the library method
                rm de_rRNA/${line}_rRNAalign.aln

                echo "=============== processing $line longest mRNA ================="
                echo "=============== processing $line =================">> b_align/mapping_report_${sort}.txt
                bowtie -p 56 -S -5 3 ${Genome}/${sort} ${line}_derRNA.fq b_align/${line}.sam 2>> b_align/mapping_report_${sort}.txt
        elif [ ${rna} = "no" ];then
                #or  don't remove rRNA
                echo "=============== processing $line longest mRNA ================="
                echo "=============== processing $line =================">> b_align/mapping_report_${sort}.txt
                bowtie -p 16 -S ${Genome}/${sort} ${sam} b_align/${line}.sam 2>> b_align/mapping_report_${sort}.txt
        fi
        cd b_align
        samtools view -bS -@ 56 ${line}.sam > ${line}.bam
        echo " ===========successfully samtobam==========="
        echo " ==================sorting=================="
        samtools sort -@ 56 ${line}.bam -o ${line}.sort.bam
        bedtools bamtobed -i ${line}.sort.bam > ${line}.bed
        sort -k1,1 -k2,2n ${line}.bed > ${line}.sort.bed
        rm ${line}.sam ${line}.bam
        cd ..
done
echo "       >>>>>>>>>>>>>>>>>>>>>>>2 longest finish<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"

## go to the fold storage .fq the genomewide processing need the result file from the step
cd "${doc}/b_align"
cp ../deadapter.txt ./
python3 ../step3_codon_pos_length.py
echo "      >>>>>>>>>>>>>>>>>>>>>>>3 step3_codon_pos_length finish<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"

#the periodicity screen
perl ../step4_trimbed.py -sp ${species} -sort ${sort_file}
cat deadapter.txt | while read line
do
bedtools genomecov -split -i ${line}_A.bed -bga -strand + -g ${Genome}/${sort}_cDNA.chrom > ${line}_A.bedgraph
done