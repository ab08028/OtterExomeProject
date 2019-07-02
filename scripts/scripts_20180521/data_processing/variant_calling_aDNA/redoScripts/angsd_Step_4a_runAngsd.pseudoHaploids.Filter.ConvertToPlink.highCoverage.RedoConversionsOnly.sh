#! /bin/bash
#$ -cwd
#$ -l h_rt=5:00:00,h_data=8G
#$ -m abe
#$ -M ab08028
#$ -N angsdStep4a
#$ -e /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -o /u/flashscratch/a/ab08028/captures/reports/angsd

######### Step 4 a: calls pseudohaploids by randomly sampling a single read at each site per individual (after applying all filters)
# based on full coverage modern + aDNA mapped to elut/mfur ######
#### run specific settings ####
trimValue=7 # set value you want to trim from either end of read (looking at mapdamage plots)
#posterior=1 # setting for angsd -doPost : 1 for using allele frequencies as prior, 2 for using a uniform prior 
#todaysdate=`date +%Y%m%d`'-highcov-pseudoHaps'
todaysdate=20190612-highcov-pseudoHaps
#### ANGSD v 0.923 ####
source /u/local/Modules/default/init/modules.sh
module load anaconda # load anaconda
source activate angsd-conda-env # activate conda env
module load plink
module load python/2.7
######### dirs and files ###########
gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptDir=$gitDir/scripts/scripts_20180521/data_processing/variant_calling_aDNA
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/aDNA-ModernComparison
bamdir=$wd/bams/
HAPdir=$wd/angsd-pseudoHaps
mkdir -p $HAPdir
mkdir -p $HAPdir/$todaysdate
outdir=$HAPdir/$todaysdate

### auxiliary scripts: 
filterHaplo=$scriptDir/filter.pseudoHaploidFile.BiallelicTransversionsOnly.py
hap2plink=/u/home/a/ab08028/klohmueldata/annabel_data/bin/angsd/misc/haploToPlink # script to convert haplo file to pseudodiploid plink file

### list of bam files to include: high coverage modern + aDNA:
elutBamList=$scriptDir/bamLists/angsd.bamList.mappedtoElutfullpaths.HighCovPlusADNAOnly.txt # list of bam files mapped to sea otter, including downsampled AND non-downsampled
mfurBamList=$scriptDir/bamLists/angsd.bamList.mappedtoMfurfullpaths.HighCovPlusADNAOnly.txt  # list of bam files mapped to ferret, including downsampled AND non-downsampled

# reference genomes:
elutRef=/u/home/a/ab08028/klohmueldata/annabel_data/sea_otter_genome/dedup_99_indexed_USETHIS/sea_otter_23May2016_bS9RH.deduped.99.fasta
mfurRef=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta

echo -e "THIS USES HIGH COVERAGE MODERN + ANCIENT ONLY\nBamLists used:\n$elutBamList\n$mfurBamList \ntrimvalue = $trimValue\nPSEUDOHAPLOIDS SAMPLING ONE RANDOM READ" > $HAPdir/$todaysdate/HIGHCOVERAGEONLY.txt


######### ANGSD settings:##############

# settings from Orlando cell paper Fages et al 2019, Cell (TAR Methods page  e14-15):
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
####### Mfur mapped bams ############
spp="mfur"
ref=$mfurRef
bamList=$mfurBamList


####### 1. get pseudo haploids ############
#angsd -nThreads 16 \
#-ref $ref \
#-bam $bamList \
#-doHaploCall 1 -doCounts 1 \
#-remove_bads 1 -uniqueOnly 1 \
#-C 50 -baq 1 -trim $trimValue -minQ 20 -minMapQ 25 \
#-out $outdir/angsdOut.mappedTo${spp}
# note skipTriallelic doesn't work in angsd here, so I do it with a custom script
# you have to do doCounts for it to work, but I'm skipping -dumpCounts 2
# don't really want the counts from this stage, it's not useful
# I want to exclude triallelic sites and filter transversions:

echo "done with angsd for $spp"
####### 2. filter so that you only have transversions (not monomorphic, or transitions) ############
# python $filterHaplo [input haplo file ] [path to output file]
# I am noting that "noRefInfo" went into calling biallelic and transversions since there could be a small number of triallelic sites slipping through
# but Fages doesn't even deal with triallelic, so I am ahead of the game already, it's fine to not exclude them with such fine resolution.
# Just note that it's just calculating biallelic and transversions on the basis of the otter sample, not based on ref genome
python $filterHaplo $outdir/angsdOut.mappedTo${spp}.haplo.gz $outdir/angsdOut.mappedTo${spp}.BiallelicTransvOnly.noRefInfo.haplo

gzip -f $outdir/angsdOut.mappedTo${spp}.BiallelicTransvOnly.noRefInfo.haplo
echo "done with filtering haplotype file"
####### 3. convert to tped format ############
# this is from the angsd git hub, not the anaconda version:
# /u/home/a/ab08028/klohmueldata/annabel_data/bin/angsd/misc/haploToPlink input.haplo.gz outputname
$hap2plink $outdir/angsdOut.mappedTo${spp}.BiallelicTransvOnly.noRefInfo.haplo.gz $outdir/angsdOut.mappedTo${spp}.BiallelicTransvOnly.noRefInfo
echo "done with converting to tped"
####### 4. convert to tped format ############
plink --tfile $outdir/angsdOut.mappedTo${spp}.BiallelicTransvOnly.noRefInfo --recode --allow-extra-chr
# and then use plink to convert tped to ped. 
echo "done with converting to ped"

####### Elut mapped bams ############
spp="elut"
ref=$elutRef
bamList=$elutBamList

####### 1. get pseudo haploids ############
#angsd -nThreads 16 \
#-ref $ref \
#-bam $bamList \
#-doHaploCall 1 -doCounts 1 \
#-remove_bads 1 -uniqueOnly 1 \
#-C 50 -baq 1 -trim $trimValue -minQ 20 -minMapQ 25 \
#-out $outdir/angsdOut.mappedTo${spp}
# note skipTriallelic doesn't work in angsd here, so I do it with a custom script
# you have to do doCounts for it to work, but I'm skipping -dumpCounts 2
# don't really want the counts from this stage, it's not useful
# I want to exclude triallelic sites and filter transversions:

echo "done with angsd for $spp"
####### 2. filter so that you only have transversions (not monomorphic, or transitions) ############
# python $filterHaplo [input haplo file ] [path to output file]
# I am noting that "noRefInfo" went into calling biallelic and transversions since there could be a small number of triallelic sites slipping through
# but Fages doesn't even deal with triallelic, so I am ahead of the game already, it's fine to not exclude them with such fine resolution.
# Just note that it's just calculating biallelic and transversions on the basis of the otter sample, not based on ref genome
python $filterHaplo $outdir/angsdOut.mappedTo${spp}.haplo.gz $outdir/angsdOut.mappedTo${spp}.BiallelicTransvOnly.noRefInfo.haplo

gzip -f $outdir/angsdOut.mappedTo${spp}.BiallelicTransvOnly.noRefInfo.haplo
echo "done with filtering haplotype file"
####### 3. convert to tped format ############
# this is from the angsd git hub, not the anaconda version:
# /u/home/a/ab08028/klohmueldata/annabel_data/bin/angsd/misc/haploToPlink input.haplo.gz outputname
$hap2plink $outdir/angsdOut.mappedTo${spp}.BiallelicTransvOnly.noRefInfo.haplo.gz $outdir/angsdOut.mappedTo${spp}.BiallelicTransvOnly.noRefInfo
echo "done with converting to tped"
####### 4. convert to tped format ############
plink --tfile $outdir/angsdOut.mappedTo${spp}.BiallelicTransvOnly.noRefInfo --recode --allow-extra-chr
# and then use plink to convert tped to ped. 
echo "done with converting to ped"


