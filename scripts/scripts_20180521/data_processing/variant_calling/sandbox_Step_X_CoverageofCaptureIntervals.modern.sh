#! /bin/bash
#$ -cwd
#$ -l h_rt=5:00:00,h_data=20G
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
module load bedtools
# header:
header=$1 # input header into file
# file locations:
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures
paleomixOutput=$wd/paleomix/${header} # specific to this header
outdir=$wd/coverage
mkdir -p $outdir
# programs
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar
# ferret reference:
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta
REFPREFIX=Mustela_putorius_furo.MusPutFur1.0.dna.toplevel

# going to be more lax on coverage requirements (10 --> 1) to account for ancient DNA.
# can always filter out later.
# note. for modern use cov = 5 to get intervals; for ancient use cov = 1. 
java -jar $GATK \
	-T DepthOfCoverage \
	-R $REFERENCE \
	-I $paleomixOutput/${header}.${REFPREFIX}.bam \
	-L $wd/captureRegions/ferret.Exon.Coordinates.0based.bed
	-o $outdir/${header}.coverage.exons
# testing this out: I think I can also do all bams at once. 
# *********** not sure if I want to do this yet ********** 
sleep 10m
