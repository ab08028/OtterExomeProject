#! /bin/bash
#$ -cwd
#$ -l h_rt=5:00:00,h_data=12G,highp
#$ -o /u/flashscratch/a/ab08028/captures/reports/fastqc
#$ -e /u/flashscratch/a/ab08028/captures/reports/fastqc
#$ -m bea
#$ -M ab08028
#$ -t 1-12

############################## ANCIENT MAPPING #######################################

i=${SGE_TASK_ID}

# program locations on Hoffman: need to use updated fastqc to deal with Novaseq
# version 0.11.7 (updated from 0.11.5)
fastqc=/u/home/a/ab08028/klohmueldata/annabel_data/bin/FastQC/fastqc

# locations (put in every script)
SCRATCH=/u/flashscratch/a/ab08028/
wd=$SCRATCH/captures
fastqs=$wd/fastqs

# go to fastqs dir
cd $fastqs

# make output directory inside the fastqs dir
mkdir -p fastqc-output


# example ancient filename: A1_Elut_CA_AN_396_SN1_S61_R1_001.fastq.gz
# ancient samples start with A...
fileR1=`ls A${i}_*_*R1*.fastq.gz` 
fileR2=`ls A${i}_*_*R2*.fastq.gz`

# header=${fileR1%_S*_R*} # save this for later

$fastqc $fileR1 -o fastqc-output
$fastqc $fileR2 -o fastqc-output


sleep 10m
