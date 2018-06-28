#! /bin/bash
#$ -cwd
#$ -l h_rt=10:00:00,h_data=21G,highp,arch=intel*
#$ -N GTgVCF
#$ -o /u/scratch2/a/ab08028/otters/reports
#$ -e /u/scratch2/a/ab08028/otters/reports
#$ -m abe
#$ -M ab08028
#$ -t 1-130

source /u/local/Modules/default/init/modules.sh
module load java/1.8.0_111
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar
# 20170808: updated reference
# 20170725: Found duplicated scaffolds in genome (redundant sequence from the assembly); consulted with Dovetail
# removed using Dedupe with 99% identity:
# $dedupe in=$ref out=./sea_otter_23May2016_bS9RH.deduped.99.fasta outd=duplicates.fa findoverlap=t cluster=t minidentity=99 csf=clusters.file
# and did the Step 4.0 to set it up as a reference genome
#REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/sea_otter_genome/dedup_99_indexed_USETHIS/sea_otter_23May2016_bS9RH.deduped.99.fasta
# 20170904: using ferret as a reference to align to. Older alignments used sea otter genome with duplicated scaffolds
# and had an issue with reference bias. Whereas this will have hets called and homs called relative to ferret

REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fa 
# have to make these before 
BAM_DIR=$1 # path to /bams dir 
IN_DIR=$BAM_DIR/ReadsFiltered
OUT_DIR=${BAM_DIR%/bams*}/gvcfs/HaplotypeCaller # instead, going to try using the gvcfs dir. 
PREFIX=$2

mkdir -p $OUT_DIR
# make interval files of each scaffold (do once)s
# index=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna_sm.toplevel.fa.fai 
#cd /u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/intervalFiles
#for i in {1..129}
#do
# the index gives the name ($1) and the length ($2)
# since bed format is [0,X) non inclusive, you can print 0 for start and length as end to get full scaffold

#awk -v i=$i '{OFS="\t";if(NR==i)print $1,0,$2}' $index > interval_${i}.bed 
#done
# then for the remaining scaffolds that are smaller, but them in one intervals file:
#tail -n +130 $index | awk '{OFS="\t";print $1,0,$2}' > interval_130.bed
# the plus character says to start printing after that line -1, so first entry should be line 100 
#
# this makes it more computationally efficient 

# make 130 files (first 129 have big scaffolds, then 130 has a big scaffold plus all the little ones)
intervalFiles=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/intervalFiles

java -jar $GATK \
-T HaplotypeCaller \
-R $REFERENCE \
-ERC BP_RESOLUTION \
-mbq 20 \
-L $intervalFiles/interval_${SGE_TASK_ID}.bed \
-out_mode EMIT_ALL_SITES \
-I $IN_DIR/${PREFIX}_Aligned_MarkDup_IndelRealigned_Filtered.bam \
-o $OUT_DIR/${PREFIX}.interval_${SGE_TASK_ID}.hapcaller.g.vcf.gz
# BP_RESOLUTION keeps all sites; do I want them in bands? how big is resulting gvcf?
# do I want to just call over my targets? 
# I decided to remove --dontUseSoftClippedBases since Jacqueline isn't sure it's necessary
