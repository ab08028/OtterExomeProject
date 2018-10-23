#! /bin/bash
#$ -cwd
#$ -l h_rt=24:00:00,h_data=8G
#$ -N blastFastqGalap
#$ -o /u/flashscratch/a/ab08028/GalapagosRails/reports
#$ -e /u/flashscratch/a/ab08028/GalapagosRails/reports
#$ -m abe
#$ -M ab08028
#$ -pe shared 10
# goal of this script: to find out what a fastq file blasts to (figure out what contamination is)

wd=/u/flashscratch/a/ab08028/GalapagosRails/
fastqs=$wd/fastqs
blastdir=/u/home/a/ab08028/bin/ncbi-blast-2.5.0+/bin/  # location of blast scripts
blastdb=/u/flashscratch/a/ab08028/NCBI_nt_db/nt # # blast db was downloaded from ncbi FTB and then all files were un-tarred on 20181022

for file in LS07_S6_L007_R1_001.fastq LS24_S7_L007_R1_001.fastq
do

input=${file%.fastq}.fasta
output=${file%.fastq}.blast.out
outdir=$wd/blast/${file%.fastq}
mkdir -p $outdir
# convert fastq --> fasta
# module load fastx_toolkit
# (unzipp fastq files ) 
# converted to fasta -- need -Q33 for illumina scores 
# fastq_to_fasta -i $fastqs/$file -o ${file%.fastq}.fasta -n -Q33


$blastdir/blastn \
-query $fastqs/$input \
-db ${blastdb} \
-outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid qlen" \
-out ${outdir}/$output \
-num_threads 10


######## I also want to blast my double capture and see what it blasts to? or my single capture? 
# blast A13 and A9 and compare? and one uncaptured SM -- A19 #

done