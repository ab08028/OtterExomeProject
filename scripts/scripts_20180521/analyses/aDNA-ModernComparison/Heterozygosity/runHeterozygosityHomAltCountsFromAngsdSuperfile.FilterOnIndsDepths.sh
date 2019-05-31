#! /bin/bash
#$ -cwd
#$ -l h_rt=40:00:00,h_data=4G,highp
#$ -m abe
#$ -M ab08028
#$ -N parseAngsdOutput
#$ -e /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -o /u/flashscratch/a/ab08028/captures/reports/angsd

####### calculate heterozygosity from posterior probabilities, exclude sites below a depth threshold on a per individual basis #######
# requires that you used -doCounts 1 -dumpCounts 2 when running ANGSD
source /u/local/Modules/default/init/modules.sh
module load python/2.7


# directories: 
gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptDir=$gitDir/scripts/scripts_20180521/
script=$scriptDir/analyses/aDNA-ModernComparison/Heterozygosity/parseBeagleSuperfile.ManyFilters.py # 20190531, this works from super file and filters on multiple things

echo $1 $2 $3 $4 $5 $6
input=$1 # path to desired superfile

sampleIDs=$2

output=$3 # path to output file
maxProbCutoff=$4 # e.g. 0.95 # this is the cutoff for the max posterior probability. If the max of one of the three GTs posteriors isn't >=
# than this cutoff, then it won't be counted for that individual. Note that it doesn't have to be the het GT that is >0.5, just one of the three
# this is to avoid cases where each of the three GTs is very close in probability, indicating low overal confidence or possibly no data
minDepthCutoff=$5 # e.g. 1 # sites that are < this threshold will not be counted toward an individuals heterozygosity
minInds=$6 # e.g. 2 # min number of individuals for a site to be worked on overall; because the prior is based on allele frequency, we don't want it to just be one individual
# note that this minInd is not just the nInd in the maf file which shows how many inds had at least 1 read. instead it's how many inds have at least minDepthCutoff reads
# which in this case is 1, but you can alter it
# new usage: 
# usage: python script.py inputFilepath sampleIDFile outputFile MaxProbCutoff PerIndividualDepthMinimum minIndsPerSite
# be careful of order, though python script will output the filters used as part of the output so you can see if anything went wrong.
python $script $input $sampleIDs $output $maxProbCutoff $minDepthCutoff $minInds

sleep 10m