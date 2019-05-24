#! /bin/bash
#$ -cwd
#$ -l h_rt=200:00:00,h_data=2G,highp
#$ -m abe
#$ -M ab08028
#$ -pe shared 16
#$ -N angsdOrlandoFullCov
#$ -e /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -o /u/flashscratch/a/ab08028/captures/reports/angsd

####### want to do full coverage modern + aDNA ######
#### ANGSD v 0.923 ####
source /u/local/Modules/default/init/modules.sh
module load anaconda # load anaconda
source activate angsd-conda-env # activate conda env
gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptDir=$gitDir/scripts/scripts_20180521/data_processing/variant_calling_aDNA
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/aDNA-ModernComparison
bamdir=$wd/bams/
GLdir=$wd/angsd-GLs
mkdir -p $GLdir
#todaysdate=`date +%Y%m%d`
todaysdate="20190523-highcov"
mkdir -p $GLdir/$todaysdate/posteriorProbabilities
outdir=$GLdir/$todaysdate/posteriorProbabilities
######### ANGSD: get posteriors ##############
# 20190523: adding -doCounts 1 -dumpCounts 2 to all; this puts out counts per individual per site (incorporating filters)
# which I then use in my heterozygosity calculations

# settings from Orlando cell paper Fages et al 2019, Cell (TAR Methods page  e14-15):
# -doMajorMinor 1 -doMaf 1 -beagleProb 1 -doPost 1 -GL 2 -minQ 20 -minMapQ 25 -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50

# doMajorMinor 1: Infer major and minor from GL
# doMaf 1: Frequency (fixed major and minor) (gets them from GLs from doMajorMinor above)
# beagleProb 1: Dump beagle style postprobs
# doPost 1: estimate the posterior genotype probability based on the allele frequency as a prior  ; ASSUMES HWE
# GL 2: 
# minQ 20 : base quality 
# mapQ 25 : mapping quality I was previously using minMapQ of 30 ; drop down to 25
# remove_bads 1 : remove bad reads
# uniqueOnly: uniquely mapping reads only
# baq 1 : lower qual scores around indels *hadn't done this previously*
# C 50: lower qual scores when there are a lot of mismatches  
### AB: I am also adding the trim 4bp on either end

# You can then get heterozygosity without the SFS ANGSD business that is buggy
# Instead, you will sum up the posterior probabilities over all sites; don't count sites with missing data in the denominator
# for each individual
# gather bams from paleomix using script gatherBamsForDownsampling.sh
# and make lists of the relevant bams: 

###### dopost1 assumes hWE


#elutBamList=$scriptDir/bamLists/angsd.bamList.mappedtoElutfullpaths.txt # list of bam files mapped to sea otter, including downsampled AND non-downsampled
#mfurBamList=$scriptDir/bamLists/angsd.bamList.mappedtoMfurfullpaths.txt  # list of bam files mapped to ferret, including downsampled AND non-downsampled
elutBamList=$scriptDir/bamLists/angsd.bamList.mappedtoElutfullpaths.HighCovPlusADNAOnly.txt # list of bam files mapped to sea otter, including downsampled AND non-downsampled
mfurBamList=$scriptDir/bamLists/angsd.bamList.mappedtoMfurfullpaths.HighCovPlusADNAOnly.txt  # list of bam files mapped to ferret, including downsampled AND non-downsampled

# references:
elutRef=/u/home/a/ab08028/klohmueldata/annabel_data/sea_otter_genome/dedup_99_indexed_USETHIS/sea_otter_23May2016_bS9RH.deduped.99.fasta
mfurRef=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta

echo "THIS USES HIGH COVERAGE MODERN + ANCIENT ONLY" > $GLdir/$todaysdate/HIGHCOVERAGEONLY.txt
# trying output in beagle format  doGlf 2
####### Mfur mapped bams ############
spp="mfur"
ref=$mfurRef
bamList=$mfurBamList

angsd -nThreads 16 \
-ref $ref \
-bam $bamList \
-GL 2 \
-doMajorMinor 1 -doMaf 1 \
-beagleProb 1 -doPost 1 \
-remove_bads 1 -uniqueOnly 1 \
-C 50 -baq 1 -trim 4 -minQ 20 -minMapQ 25 -skipTriallelic 1 \
-out $outdir/angsdOut.mappedTo${spp}.OrlandoSettings \
-doGlf 2 \
-doCounts 1 -dumpCounts 2


####### Elut mapped bams ############
spp="elut"
ref=$elutRef
bamList=$elutBamList

angsd -nThreads 16 \
-ref $ref \
-bam $bamList \
-GL 2 \
-doMajorMinor 1 -doMaf 1 \
-beagleProb 1 -doPost 1 \
-remove_bads 1 -uniqueOnly 1 \
-C 50 -baq 1 -trim 4 -minQ 20 -minMapQ 25 -skipTriallelic 1 \
-out $outdir/angsdOut.mappedTo${spp}.OrlandoSettings \
-doGlf 2 \
-doCounts 1 -dumpCounts 2



