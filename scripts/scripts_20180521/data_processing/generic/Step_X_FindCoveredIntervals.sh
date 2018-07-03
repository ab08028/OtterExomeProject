#! /bin/bash
#$ -cwd
#$ -l h_rt=20:00:00,h_data=20G,arch=intel*,highp
#$ -m bea

########### gatk find covered regions
# do this after mapping
# going to do this for each file and then can use it in variant calling.
# this will maximize called sites
# then can try to find neutral regions among these. 

# modules

source /u/local/Modules/default/init/modules.sh
module load java/1.8.0_111
module load samtools

# file locations:
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures
fastqs=$wd/fastqs
bams=$wd/bams
outdir=$wd/coveredIntervals
mkdir -p $outdir
# programs
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar
# ferret reference:
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta
# min coverage: 10 (or 5?)
# minimum map quality : 30 (JAR pipeline)
# minimum base quality: 20 (JAR pipeline)
header=$1 # input header into file
java -jar $GATK \
	-T FindCoveredIntervals \
	-R $REFERENCE \
	-I ${header}.bam \
	-cov 10 \
	-minBQ 20 \
	-minMQ 30 \
	-o $outdir/${header}.coveredIntervals.txt
# you'll then use this as -L when you call variants.