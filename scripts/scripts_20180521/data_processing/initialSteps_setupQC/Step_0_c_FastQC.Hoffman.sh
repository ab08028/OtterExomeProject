#! /bin/bash
#$ -cwd
#$ -l h_rt=5:00:00,h_data=12G
#$ -m bea

############################## ANCIENT MAPPING #######################################

header=$1

# program locations on Hoffman: need to use updated fastqc to deal with Novaseq
# version 0.11.7 (updated from 0.11.5)
fastqc=/u/home/a/ab08028/klohmueldata/annabel_data/bin/FastQC/fastqc

# locations (put in every script)
SCRATCH=/u/flashscratch/a/ab08028/
wd=$SCRATCH/captures
fastqs=$wd/fastqs
headers=$wd/samples/allElutSamples.txt # can change this to work on different samples

# make output directory
mkdir -p $wd/fastqc

# adjust if different format
fileR1=$fastqs/${header}_*_R1_*.fastq.gz
fileR2=$fastqs/${header}_*_R2_*.fastq.gz

# run fastqc
$fastqc $fileR1 -o $wd/fastqc
$fastqc $fileR2 -o $wd/fastqc

# after it's done, run multiqc separately to aggregate.
sleep 10m
