#! /bin/bash
#$ -cwd
#$ -l h_rt=5:00:00,h_data=12G
#$ -o /u/flashscratch/a/ab08028/captures/reports/step_1_fastqc
#$ -e /u/flashscratch/a/ab08028/captures/reports/step_1_fastqc
#$ -m bea
#$ -M ab08028
#$ -t 30-167

############################## ANCIENT MAPPING #######################################

i=${SGE_TASK_ID}

# program locations on Hoffman: need to use updated fastqc to deal with Novaseq
# version 0.11.7 (updated from 0.11.5)
fastqc=/u/home/a/ab08028/klohmueldata/annabel_data/bin/FastQC/fastqc

# locations (put in every script)
SCRATCH=/u/flashscratch/a/ab08028/
wd=$SCRATCH/captures
fastqs=$wd/fastqs
bams=$wd/bams

# go to fastqs dir
cd $fastqs

# make output directory
mkdir -p fastqc-output


# adjust if different format
fileR1=`ls ${i}_Elut_*R1*.fastq.gz` 
fileR2=`ls ${i}_Elut_*R2*.fastq.gz`

# run fastqc
$fastqc $fileR1 -o fastqc-output
$fastqc $fileR2 -o fastqc-output

# after it's done, run multiqc separately to aggregate.
sleep 10m
