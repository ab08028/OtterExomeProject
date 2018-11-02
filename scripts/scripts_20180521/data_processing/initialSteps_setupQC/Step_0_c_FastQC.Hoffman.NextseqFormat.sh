#! /bin/bash
#$ -cwd
#$ -l h_rt=5:00:00,h_data=12G
#$ -m bea
# next seq format for Baja samples

header=$1

# program locations on Hoffman: need to use updated fastqc to deal with Novaseq
# version 0.11.7 (updated from 0.11.5)
fastqc=/u/home/a/ab08028/klohmueldata/annabel_data/bin/FastQC/fastqc

# locations (put in every script)
SCRATCH=/u/flashscratch/a/ab08028/
wd=$SCRATCH/captures
fastqs=$wd/fastqs


# make output directory
mkdir -p $wd/fastqc

# adjust if different format

for i in {1..4}
do
for j in {1..2}
do
$fastqc $fastqs/$header_*_L00${i}_R${j}_*.fastq.gz -o $wd/fastqc
done
done
# after it's done, run multiqc separately to aggregate.
sleep 10m
