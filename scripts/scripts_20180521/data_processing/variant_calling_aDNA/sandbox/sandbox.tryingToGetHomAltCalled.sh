# Exploring ANGSD
# Trying to figure out why 1/1 GLs and GPs are so low, even for high cov aDNA
# remember that this doesnt' affect any of your other findings, so settle down. 
# maybe need to not add them up but instead 'call' them? how many have a > 0.66 frequency?
# 

 # Things to try:
 # doMajorMinor
 ######### Step 1 a-i: call GLs and GPs, counts and mafs using ANGSD based on full coverage modern + aDNA mapped to elut/mfur **using allele freqs as prior for GPS** ######
#### run specific settings ####
trimValue=7 # set value you want to trim from either end of read (looking at mapdamage plots)
posterior=1 # setting for angsd -doPost : 1 for using allele frequencies as prior, 2 for using a uniform prior 
todaysdate=`date +%Y%m%d`'-highcov-AFprior'

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


### list of bam files to include: high coverage modern + aDNA:
elutBamList=$scriptDir/bamLists/angsd.bamList.mappedtoElutfullpaths.HighCovPlusADNAOnly.txt # list of bam files mapped to sea otter, including downsampled AND non-downsampled
mfurBamList=$scriptDir/bamLists/angsd.bamList.mappedtoMfurfullpaths.HighCovPlusADNAOnly.txt  # list of bam files mapped to ferret, including downsampled AND non-downsampled

# reference genomes:
elutRef=/u/home/a/ab08028/klohmueldata/annabel_data/sea_otter_genome/dedup_99_indexed_USETHIS/sea_otter_23May2016_bS9RH.deduped.99.fasta
mfurRef=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta

echo -e "THIS USES HIGH COVERAGE MODERN + ANCIENT ONLY\nBamLists used:\n$elutBamList\n$mfurBamList \ntrimvalue = $trimValue\ndoPost posterior setting = $posterior (1 = use allele freq as prior; 2 = use uniform prior)" > $GLdir/$todaysdate/HIGHCOVERAGEONLY.txt

####### Mfur mapped bams ############
spp="mfur"
ref=$mfurRef
bamList=$mfurBamList
region="GL896898.1"
# changing:
#doMajorMinor 2 and doMaf 8 (use base counts not probs)
angsd \
-ref $ref \
-bam $bamList \
-GL 2 \
-doMajorMinor 2 -doMaf 8 \
-beagleProb 1 -doPost $posterior \
-remove_bads 1 -uniqueOnly 1 \
-C 50 -baq 1 -trim $trimValue -minQ 20 -minMapQ 25 -skipTriallelic 1 \
-out $outdir/angsdOut.mappedTo${spp}.MajorMinor2.Maf8 \
-doGlf 2 \
-doCounts 1 -dumpCounts 2 \
-r $region
