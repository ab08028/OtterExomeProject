#! /bin/bash
#$ -cwd
#$ -l h_rt=24:00:00,h_data=8G
#$ -N blastFastqOtter
#$ -o /u/flashscratch/a/ab08028/captures/reports
#$ -e /u/flashscratch/a/ab08028/captures/reports
#$ -m abe
#$ -M ab08028
#$ -pe shared 10
# goal of this script: to find out what a fastq file blasts to (figure out what contamination is)

source /u/local/Modules/default/init/modules.sh
module load fastx_toolkit

wd=/u/flashscratch/a/ab08028/captures/
fastqs=$wd/fastqs
blastdir=/u/home/a/ab08028/bin/ncbi-blast-2.5.0+/bin/  # location of blast scripts
blastdb=/u/flashscratch/a/ab08028/NCBI_nt_db/nt # # blast db was downloaded from ncbi FTB and then all files were un-tarred on 20181022

# gunzip and convert fastqs to fasta (once) - can be done in the shell 
#for file in A19_Elut_CA_SM_30_SN2_screen_S18_L007_R1_001.fastq A13_Elut_CA_AN_388_SN1_2CAP_screen_S28_L007_R1_001.fastq A9_Elut_CA_AN_388_SN1_S132_R1_001.fastq

#do
#gunzip $fastqs/$file.gz
#fastq_to_fasta -i $fastqs/$file -o $fastqs/${file%.fastq}.fasta -n -Q33
#gzip $fastqs/$file
#done

for file in A19_Elut_CA_SM_30_SN2_screen_S18_L007_R1_001.fastq A13_Elut_CA_AN_388_SN1_2CAP_screen_S28_L007_R1_001.fastq A9_Elut_CA_AN_388_SN1_S132_R1_001.fastq
do
input=${file%.fastq}.fasta
output=${file%.fastq}.blast.out
outdir=$wd/blast/${file%.fastq}
mkdir -p $outdir
 


$blastdir/blastn \
-query $fastqs/$input \
-db ${blastdb} \
-outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid qlen" \
-out ${outdir}/$output \
-num_threads 10

done