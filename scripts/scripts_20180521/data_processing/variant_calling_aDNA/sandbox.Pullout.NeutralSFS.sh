#! /bin/bash
#$ -cwd
#$ -l h_rt=10:00:00,h_data=2G,highp
#$ -m abe
#$ -M ab08028
#$ -pe shared 8
#$ -N sandbox-minDepth
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
bamdir=$wd/bams/
GLdir=$wd/angsd-GLs
SFSdir=$wd/angsd-SFS
mkdir -p $SFSdir
mkdir -p $GLdir
todaysdate=`date +%Y%m%d`
snpCutoff=1e-6
mkdir -p $GLdir/$todaysdate
mkdir -p $SFSdir/$todaysdate/perIndividual
mkdir -p $GLdir/$todaysdate/perIndividual

######## info from command line ####
bam=$1
label=$2 # ID.ref.downsample.rep#
refPrefix=$3
reference=$4

## can do realSFS in regions 
## starting with neutral (then move on to cds)
coordsDir=/u/flashscratch/a/ab08028/captures/aDNA-ModernComparison/regionCoordinates/fromModernFullDataSet/
mkdir -p $coordsDir/angsd-format
neutBed=all_8_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters.min10kb.fromExon.noCpGIsland.noRepeat.noFish.0based.sorted.merged.useThis.bed
# this is based on my full modern dataset; could probably redo based on specific coords in ancient data, but for now this is good
# convert neutBed into angsd format 
# angsd format is 1-based Chr1:40-80
# so from bed file want to do:
awk '{OFS="";print $1,":",$2+1,"-",$3}' $coordsDir/bedCoords/$neutBed > $coordsDir/angsd-format/$neutBed%.bed}.angsFmt.txt


realSFS $SFSdir/$todaysdate/perIndividual/${label}.mappedTo${refPrefix}.allSites.TransversionsOnly.saf.idx -rf $neutBed > $SFSdir/$todaysdate/perIndividual/${label}.mappedTo${refPrefix}.allSites.TransversionsOnly.saf.SFS.NEUTRAL.txt
realSFS $SFSdir/$todaysdate/perIndividual/${label}.mappedTo${refPrefix}.allSites.saf.idx $neutBed > $SFSdir/$todaysdate/perIndividual/${label}.mappedTo${refPrefix}.allSites.saf.SFS.NEUTRAL.txt

source deactivate

sleep 10m