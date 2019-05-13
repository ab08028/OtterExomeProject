#! /bin/bash
#$ -cwd
#$ -l h_rt=50:00:00,h_data=8G,highp
#$ -m abe
#$ -M ab08028
#$ -N parseBeagleHets
#$ -e /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -o /u/flashscratch/a/ab08028/captures/reports/angsd

####### calculate heterozygosity from posterior probabilities #######
source /u/local/Modules/default/init/modules.sh
module load python/2.7

angsdDate=20190511 # date you ran angsd that you're interested in 
maxProbCutoff=0.5 # this is the cutoff for the max posterior probability. If the max of one of the three GTs posteriors isn't >=
# than this cutoff, then it won't be counted for that individual. Note that it doesn't have to be the het GT that is >0.5, just one of the three
# this is to avoid cases where each of the three GTs is very close in probability, indicating low overal confidence or possibly no data


# directories: 
gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptDir=$gitDir/scripts/scripts_20180521/data_processing/variant_calling_aDNA/
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/aDNA-ModernComparison
GLdir=$wd/angsd-GLs
postDir=$GLdir/$angsdDate/posteriorProbabilities # location of your posterior probs
outdir=$wd/heterozygosityFromPosteriors/
mkdir -p $outdir

#### CAUTION CAUTION CAUTION ####
sampleIDs=$scriptDir/bamLists/SampleIDsInOrder.BeCarefulOfOrder.txt ## BE VERY CAREFUL OF THE ORDER HERE
# make sure it's IDENTICAL ORDER to the bam list you used in ANGSD, otherwise you will use the wrong individuals
# beagle files are completely order-dependent 
# can use same sample ID file for mfur and elut (bamLists were in same order )
#### CAUTION CAUTION CAUTION ####

######### mfur:
ref=mfur
input=angsdOut.mappedTo${ref}.OrlandoSettings.beagle.gprobs.gz # input file 
output=$wd/heterozygosityFromPosteriors/${input%.beagle.gprobs.gz}.hetFromPost.ProbCutoff.${maxProbCutoff}.txt
## submit parsing of beagle file:
python $script $postDir/$input $sampleIDs $output $maxProbCutoff


######### elut: 
ref=elut
input=angsdOut.mappedTo${ref}.OrlandoSettings.beagle.gprobs.gz # input file 
output=$wd/heterozygosityFromPosteriors/${input%.beagle.gprobs.gz}.hetFromPost.ProbCutoff.${maxProbCutoff}.txt
## submit parsing of beagle file:
python $script $postDir/$input $sampleIDs $output $maxProbCutoff