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
module load bedtools
# header:
#header=$1 # input header into file
bamFile=$1 # path to bam file
# file locations:
SCRATCH=/u/flashscratch/a/ab08028
#wd=$SCRATCH/captures
#paleomixOutput=$wd/paleomix/testMapping/${header}/ # for now staying in the testMapping dir; eventually going to be in a different dir! 
wd=$SCRATCH/captures/aDNA-ModernComparison
bamdir=$wd/bams/
outdir=$wd/coveredIntervals
mkdir -p $outdir
# programs
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar
# ferret reference:
#### ancient are mapped to multiple reference genomes
# mfur first: 
REFERENCE=$2 # input reference from command line to go with bamList 

# going to be more lax on coverage requirements (10 --> 1) to account for ancient DNA.
# can always filter out later.
java -jar $GATK \
	-T FindCoveredIntervals \
	-R $REFERENCE \
	-I $bamFile \
	-cov 1 \
	-minBQ 20 \
	-minMQ 30 \
	-o $outdir/${bamFile%.bam}.coveredIntervals.list

# you'll then use this as -L when you call variants.
# min coverage: 10 (or 5?) -- going down to 1 for now. 
# minimum map quality : 30 (JAR pipeline)
# minimum base quality: 20 (JAR pipeline)
# are these too stringent?

# Also want to make a bed file version
# how to: the interval list has format: 
# GL896898.1:37826-37995 (scaffold:start-stop) and is 1-based.
# So for .bed file, I need to make it 0-based on the start coord and non-inclusive for end coord, so keep the same
# so want scaffold\tstart-1\tstop for bed version (to use for other things)
########## also make a bed version: 
awk -F [:-] '{OFS="\t"; print $1,$2-1,$3}' ${bamFile%.bam}.coveredIntervals.list > ${bamFile%.bam}.coveredIntervals.list.coveredIntervals.int.bed
# want to check that nothing got made negative 1 (if it started at 0)
sed -i'' 's/-1/0/g' $outdir/${header}.${REFPREFIX}.coveredIntervals.int.bed
# update: 20180730: DO NOT MERGE. Nothing gets merged and it is just slow and sometimes errors out.

#bedtools merge -i ${bamFile%.bam}.coveredIntervals.bed > ${bamFile%.bam}.coveredIntervals.int.bed
# and want to add dots for empty six fields 
awk '{OFS="\t"; print $1,$2,$3,".",".","."}' ${bamFile%.bam}.coveredIntervals.int.bed > ${bamFile%.bam}.coveredIntervals.bed
rm ${bamFile%.bam}.coveredIntervals.int.bed

######## also want to make an ANGSD format list, which is 1-based and looks like
# chr1:1-10000 so just print $2, not $2-1
awk -F [:-] '{OFS=""; print $1,":",$2,"-",$3}' ${bamFile%.bam}.coveredIntervals.list > ${bamFile%.bam}.coveredIntervals.forANGSD.txt






sleep 10m
