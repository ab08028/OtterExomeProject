#! /bin/bash
#$ -cwd
#$ -l h_rt=5:00:00,h_data=20G,arch=intel*,highp
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

# header:
header=$1 # input header into file
# file locations:
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures
paleomixOutput=$wd/paleomix/${header} # specific to this header
outdir=$wd/coveredIntervals
mkdir -p $outdir
# programs
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar
# ferret reference:
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta
REFPREFIX=Mustela_putorius_furo.MusPutFur1.0.dna.toplevel

# header:
header=$1 # input header into file
# going to be more lax on coverage requirements (10 --> 1) to account for ancient DNA.
# can always filter out later.
java -jar $GATK \
	-T FindCoveredIntervals \
	-R $REFERENCE \
	-I $paleomixOutput/${header}.${REFPREFIX}.bam \
	-cov 1 \
	-minBQ 20 \
	-minMQ 30 \
	-o $outdir/${header}.coveredIntervals.txt
# you'll then use this as -L when you call variants.
# min coverage: 10 (or 5?) -- going down to 1 for now. 
# minimum map quality : 30 (JAR pipeline)
# minimum base quality: 20 (JAR pipeline)
# are these too stringent?
sleep 10m