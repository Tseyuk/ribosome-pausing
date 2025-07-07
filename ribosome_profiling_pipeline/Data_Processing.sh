#!/bin/sh
### created and update on 02-15-2023 by wang,zeyu as a ribosomal profiling pipeline

##mapping , likely RNA_seq analysis
#set the species
species="humanref"
doc="/disk/wzy/ribo_seq/human_ribo/gxq"
rna="yes"
position="A"
Genome="/disk/wzy/ribo_seq/human_ribo/humanref"
sort="longest"   #CDS or longest
sort_file="longest_cDNA"   #or longest_cDNA
a = "" #your adapter sequence


# go to the fold storage .fq
cd ${doc}
#find -name "*mmed.fq.gz" | xargs -P 15 -n 1 pigz -d
#ls *.fq.gz |sed 's/.fq.gz//'|uniq > deadapter.txt
mkdir b_align #save the result of mapping , longest_cDNA.fa as reference
mkdir de_rRNA

find -name "*derRNA.fq.gz" | xargs -P 15 -n 1 pigz -d

#start the mapping
cat deadapter.txt | while read line
do
        if [ ${rna} = "yes" ];then
                #remove rRNA
                echo "=============== processing $line discard rRNA ================="
                echo "=============== processing $line =================">> de_rRNA/rRNA_stats.txt
                bowtie -p 56 -v 3 -5 3 -3 6 --un=${line}_derRNA.fq ${Genome}/ribosome/rDNA ${line}_trimmed_output_29-31.fastq 2>>  de_rRNA/rRNA_stats.txt > de_rRNA/${line}_rRNAalign.aln  #note ,sometimes we need "-3 4",and another we may don't need . responding to the library method
                rm de_rRNA/${line}_rRNAalign.aln

                echo "=============== processing $line longest mRNA ================="
                echo "=============== processing $line =================">> b_align/mapping_report_${sort}.txt
                bowtie -p 56 -S -5 3 -3 6 ${Genome}/${sort} ${line}_derRNA.fq b_align/${line}.sam 2>> b_align/mapping_report_${sort}.txt
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
perl ../step3_codon_pos_length.pl
echo "      >>>>>>>>>>>>>>>>>>>>>>>3 step3_codon_pos_length finish<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"

#the periodicity screen
perl ../step4_trimbed_v2.pl -sp ${species} -sort ${sort_file} -P ${position}
cat deadapter.txt | while read line
do
bedtools genomecov -split -i ${line}_A.bed -bga -strand + -g /disk/wzy/human_ribo/humanref/${sort}_cDNA.chrom > ${line}_A.bedgraph
done

#the noscreenen
perl ../step4_trimbed_noscreen.pl -sp ${species} -p ${position} -sort ${sort_file}
echo "      >>>>>>>>>>>>>>>>>>>>>>>5 step4_trim bed finish<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"

cd "${doc}/b_align/noscreen"
cp ../deadapter.txt ./

cat deadapter.txt | while read line
do
bedtools genomecov -split -i ${line}_A.bed -bga -strand + -g /disk/wzy/human_ribo/humanref/${sort}_cDNA.chrom > ${line}_A.bedgraph
done
echo "      >>>>>>>>>>>>>>>>>>>>>>>5 convert finish<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
#the frame shift calcultion
perl ../../step6_codon_pos_123.pl -sp ${species}
echo "      >>>>>>>>>>>>>>>>>>>>>>>6 frameshift ratio finish<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
cd ../
perl ../step6_codon_pos_123.pl -sp ${species}

#maybe you can also add option width to get more accurate RPKM value (skip the specific region)
perl ../step7_cov_stat.pl -sp ${species}
Rscript ../ribo_code/step8_seq_status.R ${species} ${sort}
echo "      >>>>>>>>>>>>>>>>>>>>>>>7 calculate coverage finish<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
