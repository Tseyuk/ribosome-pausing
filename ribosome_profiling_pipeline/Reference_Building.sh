#!bin/sh
###
#creat CDS_DNA.fa and longest_cDNA.fa , and hisat2/bowtie-build index longest_cDNA.fa 
#the name are same with samples in longest_cDNA_index.R and @genome in sequence_by_refFlat.pl
name="NC10"
genome="/media/hp/disk1/wzy/Genomes"
rna="no" #set yes or no
for i in $name
do
	cd ${genome}/${name}/Sequence/WholeGenomeFasta
	## first remove all unwanted charaters from fasta file.
	sed -i -e 's/ .*//' genome.fa
	## create own refFlat.txt with gtf
	
	cd ${genome}/${name}/Genes
	## first use gtfToGenepred
	${genome}/tools/gtfToGenePred -genePredExt genes.gtf ref.txt
	#${genome}/tools/gff3ToGenePred genes.gff ref.txt
	## next move the 12th column to first
	#awk 'BEGIN{FS=OFS="\t"} {a=$12; for (i=1;i<NF; i++) $i=$(i+1); $1=a; print}' ref.txt >ref1.txt
	awk 'BEGIN {FS=OFS="\t"} {$1=$12FS$1;NF--} 1' ref.txt >ref1.txt
	## finally remove the last 4 columns

	awk 'NF{NF-=4}1' FS="\t" OFS="\t" ref1.txt > refFlat.txt
	
	##	
	cd ${genome}
	python cDNA/sequence_by_refFlat.py  -g ${name}           #note: change the species in script
	Rscript cDNA/longest_cDNA.py ${name}             #note: change the species in script
done

#creat the genome size file
cd ${genome}/${name}/Sequence/WholeGenomeFasta
samtools faidx genome.fa
cut -f1,2 genome.fa.fai > sizes.genome

#creat rRNA index
if [ ${rna}=="yes" ];then
	cd ${genome}/${name}/ribosome
	bowtie-build --threads 16 ribosomes.fa rDNA
fi
