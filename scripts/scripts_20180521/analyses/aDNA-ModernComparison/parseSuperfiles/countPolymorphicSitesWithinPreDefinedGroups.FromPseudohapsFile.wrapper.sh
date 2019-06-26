#!/bin/bash
#$ -l h_rt=24:00:00,h_data=3G
#$ -N PolymorphicCounts
#$ -cwd
#$ -m bea
#$ -o /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -e /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -M ab08028
source /u/local/Modules/default/init/modules.sh
module load python/2.7

maxMissingInds=1 # number of missing GTs allowed for a group of 3 inds (can try with 0 as well)

# directories: 
wd=$SCRATCH/captures/aDNA-ModernComparison/angsd-pseudoHaps

gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptDir=$gitDir/scripts/scripts_20180521/
script=$scriptDir/analyses/aDNA-ModernComparison/parseSuperfiles/countPolymorphicSitesWithinPreDefinedGroups.FromPseudohapsFile.py
outdir=$wd/countPolymorphicSites
mkdir -p $outdir

##### high coverage info: ######
HCDate=20190612-highcov-pseudoHaps
HCSampleIDs=$scriptDir/data_processing/variant_calling_aDNA/bamLists/SampleIDsInOrder.HighCoverageAndADNAOnly.BeCarefulOfOrder.txt
mkdir -p $outdir/$HCDate
##### high coverage info: ######
LCDate=20190612-lowcov-pseudoHaps
LCSampleIDs=$scriptDir/data_processing/variant_calling_aDNA/bamLists/SampleIDsInOrder.LowCoverageOnly.BeCarefulOfOrder.txt
mkdir -p $outdir/$LCDate

################ high coverage ##################
indir=$wd/$HCDate
for ref in mfur elut
do
# usage: python script input output sampleID maxMissing
python $script $indir/angsdOut.mappedTo${ref}.haplo.gz $outdir/$HCDate/highcov.PolymorphicCounts.mappedTo${ref}.maxMiss.${maxMissingInds}.txt $HCSampleIDs ${maxMissingInds}
# the script checks for biallelic and transversions and counts up polymorphic pseudohaploid sites
done

################ low coverage ##################
indir=$wd/$LCDate
for ref in mfur elut
do
# usage: python script input output sampleID maxMissing
python $script $indir/angsdOut.mappedTo${ref}.haplo.gz $outdir/$LCDate/lowcov.PolymorphicCounts.mappedTo${ref}.maxMiss.${maxMissingInds}.txt $LCSampleIDs ${maxMissingInds}
# the script checks for biallelic and transversions and counts up polymorphic pseudohaploid sites
done
