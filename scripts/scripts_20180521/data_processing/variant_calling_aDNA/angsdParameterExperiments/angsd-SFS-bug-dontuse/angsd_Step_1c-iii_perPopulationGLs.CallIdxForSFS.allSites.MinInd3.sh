#! /bin/bash
#$ -cwd
#$ -l h_rt=30:00:00,h_data=2G,highp
#$ -m abe
#$ -M ab08028
#$ -pe shared 16
#$ -N angsdSFS
#$ -e /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -o /u/flashscratch/a/ab08028/captures/reports/angsd
### using all sites for the SFS instead of just 1e-06 snps
# will give me a sense of what L is ... but might be some whack SNPs? #
#### ANGSD v 0.923 ####
source /u/local/Modules/default/init/modules.sh
module load anaconda # load anaconda
source activate angsd-conda-env # activate conda env
gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptDir=$gitDir/scripts/scripts_20180521/data_processing/variant_calling_aDNA/
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/aDNA-ModernComparison
bamdir=$wd/bams/
GLdir=$wd/angsd-GLs
SFSdir=$wd/angsd-SFS
mkdir -p $SFSdir
mkdir -p $GLdir
#todaysdate=`date +%Y%m%d`
todaysdate=20190506 # temporarily making date later to make sure it doesn't overwrite my c-ii results
snpCutoff=1e-6
mkdir -p $GLdir/$todaysdate
mkdir -p $SFSdir/$todaysdate
mkdir -p $GLdir/$todaysdate/perPopulation
# this is temporary -- just calling in one region to make sure angsd works
# then maybe want to call genome-wide whereever we can?
# or restrict to called regions 
testRegion="ScbS9RH_100661:10009-11075"

# gather bams from paleomix using script gatherBamsForDownsampling.sh
# and make lists of the relevant bams: 

elutBams="angsd.ancient.bamList.mappedtoElutfullpaths.txt
angsd.modernAK-downsampled.bamList.mappedtoElutfullpaths.txt
angsd.modernAK.bamList.mappedtoElutfullpaths.txt
angsd.modernCA-downsampled.bamList.mappedtoElutfullpaths.txt
angsd.modernCA.bamList.mappedtoElutfullpaths.txt" # lists of bam files mapped to sea otter, including downsampled AND non-downsampled
mfurBams="angsd.ancient.bamList.mappedtoMfurfullpaths.txt
angsd.modernAK-downsampled.bamList.mappedtoMfurfullpaths.txt
angsd.modernAK.bamList.mappedtoMfurfullpaths.txt
angsd.modernCA-downsampled.bamList.mappedtoMfurfullpaths.txt
angsd.modernCA.bamList.mappedtoMfurfullpaths.txt"  # lists of bam files mapped to ferret, including downsampled AND non-downsampled

# references:
elutRef=/u/home/a/ab08028/klohmueldata/annabel_data/sea_otter_genome/dedup_99_indexed_USETHIS/sea_otter_23May2016_bS9RH.deduped.99.fasta
mfurRef=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta

# trying output in beagle format  doGlf 2

for mfurBamList in $mfurBams
do
####### Mfur mapped bams ############
##### First want to get GLF file for each set of populations (as opposed to previous variant calling that is based on all individuals for PCA)
# this command should be EXACTLY the same as what was used in angsdVariantCalling.snpCutoff1e-6.sh
# except for: doGlf is 1 (binary glf file)
# 
# aha! problem with this nInd: 
nInd=`wc -l $scriptDir/bamLists/$mfurBamList | awk '{print $1}'`
angsd \
-GL 2 \
-trim 4 \
-nThreads 16 \
-bam $scriptDir/bamLists/$mfurBamList \
-minQ 20 -minMapQ 30 \
-skipTriallelic 1 \
-doMajorMinor 4 -ref $mfurRef \
-doGlf 1 \
-uniqueOnly 1 \
-doMaf 2 \
-out $GLdir/$todaysdate/perPopulation/${mfurBamList%.bamList.*}.mappedToMfur.allSites.minInd.$nInd \
-remove_bads 1 \
-C 50 \
-minInd $nInd 
# removed : -SNP_pval 1e-6 \

# get it without transitions: 
angsd \
-dosaf 1 \
-fold 1 \
-noTrans 1 \
-anc $mfurRef \
-ref $mfurRef \
-fai ${mfurRef}.fai \
-glf $GLdir/$todaysdate/perPopulation/${mfurBamList%.bamList.*}.mappedToMfur.allSites.minInd.$nInd.glf.gz \
-out $SFSdir/$todaysdate/${mfurBamList%.bamList.*}.mappedToMfur.allSites.minInd.$nInd.TransversionsOnly \
-nInd $nInd

# get it with transitions+transversions (better for pi estimate?): 
angsd \
-dosaf 1 \
-fold 1 \
-anc $mfurRef \
-ref $mfurRef \
-fai ${mfurRef}.fai \
-glf $GLdir/$todaysdate/perPopulation/${mfurBamList%.bamList.*}.mappedToMfur.allSites.minInd.$nInd.glf.gz \
-out $SFSdir/$todaysdate/${mfurBamList%.bamList.*}.mappedToMfur.allSites.minInd.$nInd \
-nInd $nInd


# note here nInd is 3!! not 15!!!

realSFS $SFSdir/$todaysdate/${mfurBamList%.bamList.*}.mappedToMfur.allSites.minInd.$nInd.TransversionsOnly.saf.idx > $SFSdir/$todaysdate/${mfurBamList%.bamList.*}.mappedToMfur.allSites.minInd.$nInd.TransversionsOnly.saf.SFS.txt
realSFS $SFSdir/$todaysdate/${mfurBamList%.bamList.*}.mappedToMfur.allSites.minInd.$nInd.saf.idx > $SFSdir/$todaysdate/${mfurBamList%.bamList.*}.mappedToMfur.allSites.minInd.$nInd.saf.SFS.txt

# evnetually want to use -r to make SFSes for different regions from this idx file

done


for elutBamList in $elutBams
do
####### Elut mapped bams ############
##### First want to get GLF file for each set of populations (as opposed to previous variant calling that is based on all individuals for PCA)
# this command should be EXACTLY the same as what was used in angsdVariantCalling.snpCutoff1e-6.sh
# except for: doGlf is 1 (binary glf file)


nInd=`wc -l $scriptDir/bamLists/$elutBamList | awk '{print $1}'`
angsd \
-GL 2 \
-trim 4 \
-nThreads 16 \
-bam $scriptDir/bamLists/$elutBamList \
-minQ 20 -minMapQ 30 \
-skipTriallelic 1 \
-doMajorMinor 4 -ref $elutRef \
-doGlf 1 \
-uniqueOnly 1 \
-doMaf 2 \
-out $GLdir/$todaysdate/perPopulation/${elutBamList%.bamList.*}.mappedToElut.allSites.minInd.$nInd \
-remove_bads 1 \
-C 50 \
-minInd $nInd
# removed : -SNP_pval 1e-6 \

# get it without transitions: 
angsd \
-dosaf 1 \
-fold 1 \
-noTrans 1 \
-anc $elutRef \
-ref $elutRef \
-fai ${elutRef}.fai \
-glf $GLdir/$todaysdate/perPopulation/${elutBamList%.bamList.*}.mappedToElut.allSites.minInd.$nInd.glf.gz \
-out $SFSdir/$todaysdate/${elutBamList%.bamList.*}.mappedToElut.allSites.minInd.$nInd.TransversionsOnly \
-nInd $nInd

# get it with transitions+transversions (better for pi estimate?): 
angsd \
-dosaf 1 \
-fold 1 \
-anc $elutRef \
-ref $elutRef \
-fai ${elutRef}.fai \
-glf $GLdir/$todaysdate/perPopulation/${elutBamList%.bamList.*}.mappedToElut.minInd.$nInd.allSites.glf.gz \
-out $SFSdir/$todaysdate/${elutBamList%.bamList.*}.mappedToElut.minInd.$nInd.allSites \
-nInd $nInd


# note here nInd is 3!! not 15!!!

realSFS $SFSdir/$todaysdate/${elutBamList%.bamList.*}.mappedToElut.allSites.minInd.$nInd.TransversionsOnly.saf.idx > $SFSdir/$todaysdate/${elutBamList%.bamList.*}.mappedToElut.allSites.minInd.$nInd.TransversionsOnly.saf.SFS.txt
realSFS $SFSdir/$todaysdate/${elutBamList%.bamList.*}.mappedToElut.allSites.minInd.$nInd.saf.idx > $SFSdir/$todaysdate/${elutBamList%.bamList.*}.mappedToElut.allSites.minInd.$nInd.saf.SFS.txt

# evnetually want to use -r to make SFSes for different regions from this idx file
done

# can make neutral/cds sfses using the -r regions
# need to get my mfur regions in angsd format (hold off for now) 
# realSFS $GLdir/$todaysdate/angsdOut.mappedToMfur.${snpCutoff}.snpsOnly.saf.idx -r $regions

# 20190502 -- was run without doDepth or doCount
# 201090503 -- going to add more things: doDepth/doCount to get depth per sample
# Adding more filtering :
# remove_bads and -C50

source deactivate

sleep 10m