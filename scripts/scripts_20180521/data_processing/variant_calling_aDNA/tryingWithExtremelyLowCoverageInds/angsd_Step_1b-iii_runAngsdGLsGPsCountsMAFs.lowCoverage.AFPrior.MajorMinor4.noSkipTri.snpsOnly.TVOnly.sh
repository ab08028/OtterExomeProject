#! /bin/bash
#$ -cwd
#$ -l h_rt=200:00:00,h_data=2G,highp
#$ -m abe
#$ -M ab08028
#$ -pe shared 16
#$ -N angsdStep1biii
#$ -e /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -o /u/flashscratch/a/ab08028/captures/reports/angsd
############# snps only; transversions only  ############
######### Step 1 b-i: call GLs and GPs, counts and mafs using ANGSD based on downsampled modern + aDNA mapped to elut/mfur **using allele freqs as prior for GPS** ######
#### run specific settings ####
trimValue=7 # set value you want to trim from either end of read (looking at mapdamage plots)
posterior=1 # setting for angsd -doPost : 1 for using allele frequencies as prior, 2 for using a uniform prior 
snpCutoff=1e-06
#todaysdate=`date +%Y%m%d`'-lowcov-AFprior'

todaysdate='20191201-lowcov-AFprior-MajorMinor4-IncludesMoreInds-notAllIndelRealignedYet'
#### ANGSD v 0.923 ####
source /u/local/Modules/default/init/modules.sh
module load anaconda # load anaconda
source activate angsd-conda-env # activate conda env

######### dirs and files ###########
gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptDir=$gitDir/scripts/scripts_20180521/data_processing/variant_calling_aDNA
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/aDNA-ModernComparison
bamdir=$wd/bams/
GLdir=$wd/angsd-GLs
mkdir -p $GLdir
mkdir -p $GLdir/$todaysdate
outdir=$GLdir/$todaysdate

### list of bam files to include: downsampled coverage modern + aDNA:
elutBamList=$scriptDir/bamLists/angsd.bamList.LowCoverageOnly.mappedtoElutfullpaths.PLUSEXTRAADNA.notIndelRealYet.txt

# reference genomes:
elutRef=/u/home/a/ab08028/klohmueldata/annabel_data/sea_otter_genome/dedup_99_indexed_USETHIS/sea_otter_23May2016_bS9RH.deduped.99.fasta

echo -e "THIS USES LOW COVERAGE DOWNSAMPLED MODERN + ANCIENT ONLY\nBamLists used:\n$elutBamList\n$mfurBamList \ntrimvalue = $trimValue\ndoPost posterior setting = $posterior (1 = use allele freq as prior; 2 = use uniform prior)" > $GLdir/$todaysdate/LOWCOVERAGEONLY.txt

######### ANGSD settings:##############

# settings from Orlando cell paper Fages et al 2019, Cell (TAR Methods page  e14-15):
# ****** 20190701 fix: use doMajorMinor 4 instead of 1, otherwise you miss the hom-alt calls! ****** # 
# *** also not skipping triallelic, will do that myself as needed *** # 
# -doMajorMinor 1 -doMaf 1 -beagleProb 1 -doPost 1 -GL 2 -minQ 20 -minMapQ 25 -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50
# doMajorMinor 1: Infer major and minor from GL
# doMaf 1: Frequency (fixed major and minor) (gets them from GLs from doMajorMinor above)
# beagleProb 1: Dump beagle style postprobs
# doPost 1: estimate the posterior genotype probability based on the allele frequency as a prior  ; NOTE IT ASSUMES HWE
# GL 2: 
# minQ 20 : base quality 
# mapQ 25 : mapping quality I was previously using minMapQ of 30 ; drop down to 25
# remove_bads 1 : remove bad reads
# uniqueOnly: uniquely mapping reads only
# baq 1 : lower qual scores around indels *hadn't done this previously*
# C 50: lower qual scores when there are a lot of mismatches  

### AB: I am also adding the trim Xbp on either end; previously was trimming 4bp, based on mapdamage I want to do 7bp 
# I'm also adding -doCounts 1 -dumpCounts 2 to all; this puts out counts per individual per site (incorporating filters)
# which I then use in my heterozygosity calculations

###### NOTE : doSAF and realSFS are extremely buggy in ANGSD-- don't work with scaffold/capture data when some scaffold are missing data
# So you don't want to use those; use my custom downstream scripts instead.

# trying output in beagle format  doGlf 2

####### Elut mapped bams ############
spp="elut"
ref=$elutRef
bamList=$elutBamList

angsd -nThreads 16 \
-ref $ref \
-bam $bamList \
-GL 2 \
-doMajorMinor 4 -doMaf 1 \
-beagleProb 1 -doPost $posterior \
-remove_bads 1 -uniqueOnly 1 \
-C 50 -baq 1 -trim $trimValue -minQ 20 -minMapQ 25 \
-out $outdir/angsdOut.mappedTo${spp}.${snpCutoff}.snpsOnly.TransvOnly \
-doGlf 2 \
-doCounts 1 -dumpCounts 2 \
-SNP_pval $snpCutoff \
-rmTrans 1

source deactivate


