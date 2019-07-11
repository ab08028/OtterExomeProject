#### wrapper to submit:
# directories: 
#### dirs ####
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/aDNA-ModernComparison
GLdir=$wd/angsd-GLs
gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptDir=$gitDir/scripts/scripts_20180521/
script=$scriptDir/analyses/aDNA-ModernComparison/Heterozygosity/runHeterozygosityHomAltCountsFromAngsdSuperfile.FilterOnIndsDepths.sh
# three filters:
# I like this more stringent maxProbCutoff filter; it doesn't affect that many sites for GPs;
maxProbCutoff=0.95 # this is the cutoff for the max posterior probability. If the max of one of the three GTs posteriors isn't >=
#maxProbCutoff=0.5 # also trying 0.5
# than this cutoff, then it won't be counted for that individual. Note that it doesn't have to be the het GT that is >0.5, just one of the three
# this is to avoid cases where each of the three GTs is very close in probability, indicating low overal confidence or possibly no data
minDepthCutoff=1 # sites that are < this threshold will not be counted toward an individuals heterozygosity
# want to try a range minInds
#minIndsOptions="2 3 5 9" # options for minInds ; can try multiple combinations to see how you lose sites
minIndsOptions=2 # min number of individuals for a site to be worked on overall; because the prior is based on allele frequency, we don't want it to just be one individual
# note that this minInd is not just the nInd in the maf file which shows how many inds had at least 1 read. instead it's how many inds have at least minDepthCutoff reads
# which in this case is 1, but you can alter it (sites with < minInds will not be counted for ANY individuals even if they have data.)
hcDates="20190701-highcov-AFprior-MajorMinor4" # high cov dates ; skipping UNIF prior because it is garbage: 20190524-highcov-UNIFprior
lcDates="20190701-lowcov-AFprior-MajorMinor4"  # low cov dates
# high coverage:
# for ref in $refs
refs="mfur elut" 

################################ GPs #########################
type="GPs"
for ref in $refs
do
superfile=angsdOut.mappedTo${ref}.superfile.${type}.mafs.counts.0based.bed.gz # name of superfile; can be GLs or GPs; can be whole genome, cds, neutral, etc. as long as it's in superfile format (bed fmt + maf file + GP or GL + counts)
######### high coverage: #########
for angsdDate in $hcDates
do
for minInds in $minIndsOptions
do
# these are the high coverage sample IDs:
sampleIDs=$scriptDir/data_processing/variant_calling_aDNA/bamLists/SampleIDsInOrder.HighCoverageAndADNAOnly.BeCarefulOfOrder.txt

# get dirs:
indir=$GLdir/$angsdDate 
outdir=$wd/heterozygosityFromPosteriors/$angsdDate
mkdir -p $outdir
output=$outdir/${superfile%.mafs.counts.0based.bed.gz}.hetHomTotals.ProbCutoff.${maxProbCutoff}.DepthCutoff.${minDepthCutoff}.minInd.${minInds}.${angsdDate}.txt
# order is
# qsub -N name script -v "inputsuperfile sampleIDFile outputFile maxProbCutoff minDepthCutoff minIndCutoff" # use -F "" to pass variables to pass into runHeterozygosityHomAltCountsFromAngsdSuperfile.FilterOnIndsDepths.sh
qsub -N parseHC${ref}${type}${minInds} $script $indir/$superfile $sampleIDs $output $maxProbCutoff $minDepthCutoff $minInds
done
done
# -N parseHC${ref}${type}${minInds} 
######### low coverage ##############
for angsdDate in $lcDates
do
for minInds in $minIndsOptions
do
# these are the downsampled low coverage sample IDs:
sampleIDs=$scriptDir/data_processing/variant_calling_aDNA/bamLists/SampleIDsInOrder.LowCoverageOnly.BeCarefulOfOrder.txt


# get dirs:
indir=$GLdir/$angsdDate 
outdir=$wd/heterozygosityFromPosteriors/$angsdDate
mkdir -p $outdir
output=$outdir/${superfile%.mafs.counts.0based.bed.gz}.hetHomTotals.ProbCutoff.${maxProbCutoff}.DepthCutoff.${minDepthCutoff}.minInd.${minInds}.${angsdDate}.txt


qsub -N parseLC${ref}${type}${minInds} $script $indir/$superfile $sampleIDs $output $maxProbCutoff $minDepthCutoff $minInds
done

done
done

############################## GLs: may want different settings ########################################
######### skipping GLs for now because may want different filter settings #######################
