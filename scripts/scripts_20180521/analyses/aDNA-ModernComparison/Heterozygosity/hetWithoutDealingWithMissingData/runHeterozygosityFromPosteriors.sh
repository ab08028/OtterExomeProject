#! /bin/bash
#$ -cwd
#$ -l h_rt=100:00:00,h_data=8G,highp
#$ -m abe
#$ -M ab08028
#$ -N hetPerIndMissingData
#$ -e /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -o /u/flashscratch/a/ab08028/captures/reports/angsd

####### calculate heterozygosity from posterior probabilities #######
source /u/local/Modules/default/init/modules.sh
module load python/2.7

maxProbCutoff=0.5 # this is the cutoff for the max posterior probability. If the max of one of the three GTs posteriors isn't >=
# than this cutoff, then it won't be counted for that individual. Note that it doesn't have to be the het GT that is >0.5, just one of the three
# this is to avoid cases where each of the three GTs is very close in probability, indicating low overal confidence or possibly no data


# directories: 
gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptDir=$gitDir/scripts/scripts_20180521/
script=$scriptDir/analyses/aDNA-ModernComparison/Heterozygosity/hetWithoutDealingWithMissingData/parseBeaglePosteriors.py
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/aDNA-ModernComparison
GLdir=$wd/angsd-GLs

############# all samples #######################
# angsdDate=20190511 # date you ran angsd that you're interested in 
# postDir=$GLdir/$angsdDate/posteriorProbabilities # location of your posterior probs
# outdir=$wd/heterozygosityFromPosteriors/$angsdDate
# mkdir -p $outdir
# #### CAUTION CAUTION CAUTION ####
# sampleIDs=$scriptDir/data_processing/variant_calling_aDNA/bamLists/SampleIDsInOrder.BeCarefulOfOrder.txt ## BE VERY CAREFUL OF THE ORDER HERE
# # make sure it's IDENTICAL ORDER to the bam list you used in ANGSD, otherwise you will use the wrong individuals
# # beagle files are completely order-dependent 
# # can use same sample ID file for mfur and elut (bamLists were in same order )
# #### CAUTION CAUTION CAUTION ####
# 
# ######### mfur:
# ref=mfur
# input=angsdOut.mappedTo${ref}.OrlandoSettings.beagle.gprobs.gz # input file 
# output=$outdir/${input%.beagle.gprobs.gz}.hetFromPost.ProbCutoff.${maxProbCutoff}.${angsdDate}.txt
# ## submit parsing of beagle file:
# python $script $postDir/$input $sampleIDs $output $maxProbCutoff
# 
# 
# ######### elut: 
# ref=elut
# input=angsdOut.mappedTo${ref}.OrlandoSettings.beagle.gprobs.gz # input file 
# output=$outdir/${input%.beagle.gprobs.gz}.hetFromPost.ProbCutoff.${maxProbCutoff}.${angsdDate}.txt
# ## submit parsing of beagle file:
# python $script $postDir/$input $sampleIDs $output $maxProbCutoff


############# high coverage +aDNA only (with and wihtout minind 5) #######################
sampleIDs=$scriptDir/data_processing/variant_calling_aDNA/bamLists/SampleIDsInOrder.HighCoverageAndADNAOnly.BeCarefulOfOrder.txt

#for angsdDate in 20190513-highcov-minInd 20190513-highcov
# try the neut-only (first run, didn't have the Counts file)
for angsdDate in 20190521-highcov-neutOnly
do
postDir=$GLdir/$angsdDate/posteriorProbabilities # location of your posterior probs
outdir=$wd/heterozygosityFromPosteriors/$angsdDate
mkdir -p $outdir

ref=mfur
input=angsdOut.mappedTo${ref}.OrlandoSettings.beagle.gprobs.gz # input file 
output=$outdir/${input%.beagle.gprobs.gz}.hetFromPost.ProbCutoff.${maxProbCutoff}.${angsdDate}.txt
python $script $postDir/$input $sampleIDs $output $maxProbCutoff

#ref=elut
#input=angsdOut.mappedTo${ref}.OrlandoSettings.beagle.gprobs.gz # input file 
#output=$outdir/${input%.beagle.gprobs.gz}.hetFromPost.ProbCutoff.${maxProbCutoff}.${angsdDate}.txt
#python $script $postDir/$input $sampleIDs $output $maxProbCutoff
done


########## low coverage only ###########
sampleIDs=$scriptDir/data_processing/variant_calling_aDNA/bamLists/SampleIDsInOrder.LowCoverageOnly.BeCarefulOfOrder.txt


#for angsdDate in 20190513-lowcov-minInd 20190513-lowcov
for angsdDate in 20190521-lowcov-neutOnly
do
postDir=$GLdir/$angsdDate/posteriorProbabilities # location of your posterior probs
outdir=$wd/heterozygosityFromPosteriors/$angsdDate
mkdir -p $outdir

ref=mfur
input=angsdOut.mappedTo${ref}.OrlandoSettings.beagle.gprobs.gz # input file 
output=$outdir/${input%.beagle.gprobs.gz}.hetFromPost.ProbCutoff.${maxProbCutoff}.${angsdDate}.txt
python $script $postDir/$input $sampleIDs $output $maxProbCutoff

#ref=elut
#input=angsdOut.mappedTo${ref}.OrlandoSettings.beagle.gprobs.gz # input file 
#output=$outdir/${input%.beagle.gprobs.gz}.hetFromPost.ProbCutoff.${maxProbCutoff}.${angsdDate}.txt
#python $script $postDir/$input $sampleIDs $output $maxProbCutoff
done
