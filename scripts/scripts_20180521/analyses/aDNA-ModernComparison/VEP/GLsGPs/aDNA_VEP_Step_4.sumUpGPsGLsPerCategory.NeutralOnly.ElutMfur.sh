#! /bin/bash
#$ -cwd
#$ -l h_rt=23:00:00,h_data=6G
#$ -m abe
#$ -M ab08028
#$ -N getHomAltSumsPerCat
#$ -e /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -o /u/flashscratch/a/ab08028/captures/reports/angsd


######## goal of this script is to sum up derived GPs in neutral regions (won't be the same neutral regions)
# for elut and mfur
# slight confounding issue: elut will be direct target regions, whereas mfur may have some off target stuff - so that might bias things
# hope is that they don't show a strong reference bias in derived allele fractions once corrected for number of sites
# however, a couple things could be affecting it
# 1. off-target regions -- the mfur neutral regions are determined independently of the elut regions just based on where there was coverage
# 2. biased drop out -- the capture targets that don't map as well to mfur may be more highly diverged, so lower cov individuals show fewer derived sites

# hopefully this is just a depth issue rather than a depthXref issue, but we shall see.
gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptDir=$gitDir/scripts/scripts_20180521/
script=$scriptDir/analyses/aDNA-ModernComparison/Heterozygosity/parseBeagleSuperfile.CalculateHeterozygosity.ManyFilters.py
GLdir=/u/flashscratch/a/ab08028/captures/aDNA-ModernComparison/angsd-GLs/
outdir=/u/flashscratch/a/ab08028/captures/aDNA-ModernComparison/VEP/sumGPsGLsNeutralSitesOnly
mkdir -p $outdir 
# moving this inside the VEP dir
hcdates="20190701-highcov-AFprior-MajorMinor4"
lcdates="20190701-lowcov-AFprior-MajorMinor4"
refs="elut mfur"
type=GPs # using GPs for now
categories="synonymous missense stopgained"
maxProbCutoff=0.95 # e.g. 0.95 # this is the cutoff for the max posterior probability. If the max of one of the three GTs posteriors isn't >=
# than this cutoff, then it won't be counted for that individual. Note that it doesn't have to be the het GT that is >0.5, just one of the three
# this is to avoid cases where each of the three GTs is very close in probability, indicating low overal confidence or possibly no data
minDepthCutoffs="1 2 4" # minimum only 1 read (maybe try raising this)
minInds=1 # okay to just have one individual with sequence

################################ high coverage ###############################
sampleIDs=$scriptDir/data_processing/variant_calling_aDNA/bamLists/SampleIDsInOrder.HighCoverageAndADNAOnly.BeCarefulOfOrder.txt # high cov
for ref in refs
do
basename=angsdOut.mappedTo${ref}
for minDepthCutoff in $minDepthCutoffs
do
for angsdDate in $hcdates
do
indir=$GLdir/$angsdDate 
mkdir -p $outdir/$angsdDate
neutralSites=${basename}.superfile.${type}.mafs.counts.neutralOnly.0based.bed.gz
output=${basename}.hetHomTotals.${type}.ProbCutoff.${maxProbCutoff}.DepthCutoff.${minDepthCutoff}.minInd.${minInds}.${angsdDate}.NeutralOnly.txt
# get total callable cds sites with same filters: 
python $script $indir/$neutralSites $sampleIDs $outdir/$angsdDate/$output $maxProbCutoff $minDepthCutoff $minInds

done
done
done

############################ low coverage #########################################  
sampleIDs=$scriptDir/data_processing/variant_calling_aDNA/bamLists/SampleIDsInOrder.LowCoverageOnly.BeCarefulOfOrder.txt # low cov
for ref in refs
do
basename=angsdOut.mappedTo${ref}
for minDepthCutoff in $minDepthCutoffs
do
for angsdDate in $lcdates
do
indir=$GLdir/$angsdDate 
mkdir -p $outdir/$angsdDate
neutralSites=${basename}.superfile.${type}.mafs.counts.neutralOnly.0based.bed.gz
output=${basename}.hetHomTotals.${type}.ProbCutoff.${maxProbCutoff}.DepthCutoff.${minDepthCutoff}.minInd.${minInds}.${angsdDate}.NeutralOnly.txt
# get total callable cds sites with same filters: 
python $script $indir/$neutralSites $sampleIDs $outdir/$angsdDate/$output $maxProbCutoff $minDepthCutoff $minInds

done
done
done