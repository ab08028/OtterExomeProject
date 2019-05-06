#! /bin/bash
#$ -cwd
#$ -l h_rt=10:00:00,h_data=2G,highp
#$ -m abe
#$ -M ab08028
#$ -pe shared 16
#$ -N angsdDoDepth
#$ -e /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -o /u/flashscratch/a/ab08028/captures/reports/angsd

#### ANGSD v 0.923 ####
source /u/local/Modules/default/init/modules.sh
module load anaconda # load anaconda
source activate angsd-conda-env # activate conda env
gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptDir=$gitDir/scripts/scripts_20180521/data_processing/variant_calling_aDNA/
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/aDNA-ModernComparison
#bamdir=$wd/bams/
outdir=$wd/angsd-Depths
mkdir -p $outdir
#GLdir=$wd/angsd-GLs
#mkdir -p $GLdir
# todaysdate=20190503
#todaysdate=`date +%Y%m%d`
#mkdir -p $GLdir/$todaysdate
# this is temporary -- just calling in one region to make sure angsd works
# then maybe want to call genome-wide whereever we can?
# or restrict to called regions 
testRegion="ScbS9RH_100661:10009-11075"

# gather bams from paleomix using script gatherBamsForDownsampling.sh
# and make lists of the relevant bams: 

elutBamList=$scriptDir/bamLists/angsd.bamList.mappedtoElutfullpaths.txt # list of bam files mapped to sea otter, including downsampled AND non-downsampled
mfurBamList=$scriptDir/bamLists/angsd.bamList.mappedtoMfurfullpaths.txt  # list of bam files mapped to ferret, including downsampled AND non-downsampled

# references:
elutRef=/u/home/a/ab08028/klohmueldata/annabel_data/sea_otter_genome/dedup_99_indexed_USETHIS/sea_otter_23May2016_bS9RH.deduped.99.fasta
mfurRef=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta

# trying output in beagle format  doGlf 2
####### Mfur mapped bams ############

angsd \
-trim 4 \
-nThreads 16 \
-bam $mfurBamList \
-minQ 20 -minMapQ 30 \
-uniqueOnly 1 \
-out $outdir/angsd.Depths.mappedToMFur \
-remove_bads 1 \
-C 50 \
-doDepth 1 -doCounts 1 \
-nInd 15 \
-ref $mfurRef

####### Elut mapped bams ############

angsd \
-trim 4 \
-nThreads 16 \
-bam $elutBamList \
-minQ 20 -minMapQ 30 \
-uniqueOnly 1 \
-out $outdir/angsd.Depths.mappedToElut \
-remove_bads 1 \
-C 50 \
-doDepth 1 -doCounts 1 \
-nInd 15 \
-ref $elutRef

source deactivate
