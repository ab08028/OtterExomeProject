#! /bin/bash
#$ -cwd
#$ -l h_rt=23:00:00,h_data=6G
#$ -m abe
#$ -M ab08028
#$ -N getHomAltSumsPerCat
#$ -e /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -o /u/flashscratch/a/ab08028/captures/reports/angsd


######### run het script on VEP output #########

# want to total up a few things with the same filters:
# total 'callable' cds sites from angsdOut.mappedTomfur.superfile.GPs.mafs.counts.cdsOnly.0based.bed.gz
# and then sums of GPs for 1/1 and 0/1 sites for syn, mis, and sg
gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptDir=$gitDir/scripts/scripts_20180521/
script=$scriptDir/analyses/aDNA-ModernComparison/Heterozygosity/parseBeagleSuperfile.CalculateHeterozygosity.ManyFilters.py
GLdir=/u/flashscratch/a/ab08028/captures/aDNA-ModernComparison/angsd-GLs/
outdir=/u/flashscratch/a/ab08028/captures/aDNA-ModernComparison/VEP/sumGPsGLsPerVEPCategory # moving this inside the VEP dir
hcdates="20190701-highcov-AFprior-MajorMinor4"
lcdates="20190701-lowcov-AFprior-MajorMinor4"
ref=mfur #only mfur for vep stuff
type=GPs # using GPs for now
basename=angsdOut.mappedTo${ref}
categories="synonymous missense stopgained"


maxProbCutoff=0.95 # e.g. 0.95 # this is the cutoff for the max posterior probability. If the max of one of the three GTs posteriors isn't >=
# than this cutoff, then it won't be counted for that individual. Note that it doesn't have to be the het GT that is >0.5, just one of the three
# this is to avoid cases where each of the three GTs is very close in probability, indicating low overal confidence or possibly no data
minDepthCutoffs="1 2 4" # minimum only 1 read (maybe try raising this)
minInds=1 # okay to just have one individual with sequence
################################ high coverage ###############################
sampleIDs=$scriptDir/data_processing/variant_calling_aDNA/bamLists/SampleIDsInOrder.HighCoverageAndADNAOnly.BeCarefulOfOrder.txt # high cov
for minDepthCutoff in $minDepthCutoffs
do
for angsdDate in $hcdates
do
indir=$GLdir/$date 
mkdir -p $outdir/$date
cds=${basename}.superfile.${type}.mafs.counts.cdsOnly.0based.bed.gz
output=${basename}.hetHomTotals.ProbCutoff.${maxProbCutoff}.DepthCutoff.${minDepthCutoff}.minInd.${minInds}.${angsdDate}.CDS.txt
# get total callable cds sites with same filters: 
python $script $indir/$cds $sampleIDs $outdir/$date/$output $maxProbCutoff $minDepthCutoff $minInds

for category in $categories
do
input=cdsPerCategoryFromVEP/${basename}.superfile.${type}.mafs.counts.0based.${category}.bed.gz
output=${basename}.hetHomTotals.ProbCutoff.${maxProbCutoff}.DepthCutoff.${minDepthCutoff}.minInd.${minInds}.${angsdDate}.${category}.txt

python $script $indir/$input $sampleIDs $outdir/$date/$output $maxProbCutoff $minDepthCutoff $minInds

done
done
done


############################ low coverage #########################################  
sampleIDs=$scriptDir/data_processing/variant_calling_aDNA/bamLists/SampleIDsInOrder.LowCoverageOnly.BeCarefulOfOrder.txt # low cov
for minDepthCutoff in $minDepthCutoffs
do
for angsdDate in $lcdates
do
indir=$GLdir/$date 
mkdir -p $outdir/$date
cds=${basename}.superfile.${type}.mafs.counts.cdsOnly.0based.bed.gz
output=${basename}.hetHomTotals.ProbCutoff.${maxProbCutoff}.DepthCutoff.${minDepthCutoff}.minInd.${minInds}.${angsdDate}.CDS.txt
# get total callable cds sites with same filters: 
python $script $indir/$cds $sampleIDs $outdir/$date/$output $maxProbCutoff $minDepthCutoff $minInds

for category in $categories
do
input=cdsPerCategoryFromVEP/${basename}.superfile.${type}.mafs.counts.0based.${category}.bed.gz
output=${basename}.hetHomTotals.ProbCutoff.${maxProbCutoff}.DepthCutoff.${minDepthCutoff}.minInd.${minInds}.${angsdDate}.${category}.txt

python $script $indir/$input $sampleIDs $outdir/$date/$output $maxProbCutoff $minDepthCutoff $minInds

done
done
done
# want to gzip beds afterward -- or does script only work on gzipped files?