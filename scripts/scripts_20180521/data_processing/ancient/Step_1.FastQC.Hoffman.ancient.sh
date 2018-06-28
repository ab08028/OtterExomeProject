#! /bin/bash
#$ -cwd
#$ -l h_rt=5:00:00,h_data=12G,highp
#$ -o /u/flashscratch/a/ab08028/captures/reports
#$ -e /u/flashscratch/a/ab08028/captures/reports
#$ -m bea
#$ -M ab08028
#$ -t 1-12

############################## ANCIENT MAPPING #######################################
###### Step 0.a : On Sirius Or Hoffman: fastQC
# Run from fastqs directory 
### can do this on Hoffman because not generating new fastqcs after TrimGalore

i=${SGE_TASK_ID}

# program locations on Hoffman:
fastqc=/u/home/a/ab08028/klohmueldata/annabel_data/bin/FastQC/fastqc

# locations (put in every script)
SCRATCH=/u/flashscratch/a/ab08028/
wd=$SCRATCH/captures
fastqs=$wd/fastqs
bams=$wd/bams

cd $fastqs

# make output directory
mkdir -p fastqc-output


# example ancient filename: A1_Elut_CA_AN_396_SN1_S61_R1_001.fastq.gz
# ancient samples start with A.
fileR1=`ls A${i}_Elut_*R1*.fastq.gz` 
fileR2=`ls A${i}_Elut_*R2*.fastq.gz`

# header=${fileR1%_S*_R*} # save this for later

$fastqc $fileR1 -o fastqc-output
$fastqc $fileR2 -o fastqc-output


sleep 10m
