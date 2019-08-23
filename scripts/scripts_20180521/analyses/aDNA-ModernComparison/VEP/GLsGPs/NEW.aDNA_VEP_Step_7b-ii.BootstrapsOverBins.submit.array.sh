#! /bin/bash
#$ -cwd
#$ -l h_rt=100:00:00,h_data=30G,highp
#$ -m abe
#$ -M ab08028
#$ -N bootstrapVEPResults
#$ -e /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -o /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -V
source /u/local/Modules/default/init/modules.sh
module load R/3.5.1  


gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptDir=$gitDir/scripts/scripts_20180521/analyses/aDNA-ModernComparison/VEP/GLsGPs/
script=NEW.aDNA_VEP_Step_7b-ii-alt.bootstrapBinsPerGroup.R

# this script will get point estimates based on the bins passing filters and generate bootstraps

numBoots=1000

dates="20190701-lowcov-AFprior-MajorMinor4 20190701-highcov-AFprior-MajorMinor4"
# just need these min variables to pull the correct avg site counts file:
minDepth=2 # I calculated totals with 1, 2 and 4. 2 lowers values of homAlt compared to 1. I think 2 is fitting since 1 read might seem homAlt but not be. 
minGP=0.95
minInd=1

SCRATCH=/u/flashscratch/a/ab08028/
avgSitesFile=$SCRATCH/captures/aDNA-ModernComparison/VEP/sumGPsGLsPerVEPCategory/AVERAGECALLEDSITES.allInds.HighCov.LowCov.minDepth.${minDepth}.minInd.${minInd}.minGP.${minGP}.txt
# contains both dates
ref="mfur"
basename=angsdOut.mappedTo${ref}
type="GPs" # for now

# loop through dates:
for angsdDate in $dates
do
echo $angsdDate
outPREFIX=${basename}.Bootstraps.${type}.ProbCutoff.${minGP}.DepthCutoff.${minDepth}.minInd.${minInd}.${angsdDate}
indir=$SCRATCH/captures/aDNA-ModernComparison/VEP/compareMisSynDists_withBootstraps/$angsdDate
#infile=$indir/testingScript.CDS.bed.gz
infile=$indir/${basename}.Bins.${type}.ProbCutoff.${minGP}.DepthCutoff.${minDepth}.minInd.${minInd}.${angsdDate}.Ind.${indNum}.sumsPerBin.txt # result of step 7b-i 
 # need to fix header issue
outdir=$indir/PointEstsPlusBootstraps
mkdir -p $outdir
# file with avg sites to draw (depends on minDepth,minGP, minInd) and is a result of the script: /Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/analyses/aDNA-ModernComparison/VEP/GLsGPs/NEW_aDNA_VEP_Step_7b-i-alt.setUpBootstrapBinsPerAncientModern.submit.array.sh
# get count for this date:
# parameters for R script:
# option_list = list(
#   make_option(c("--infile"), type="character", default=NULL, 
#               help="path to your input file result of step 7b-i which generates per indivdiual bins across the genome with the GPs of each category (syn, mis, sg) summed up. This is the df indAllBins from previous step", metavar="file"),
#   make_option(c("--outdir"), type="character", default=NULL, 
#               help="path to outdir", metavar="dir"),
#   make_option(c("--outPREFIX"), type="character", default=NULL, 
#               help="outfilePrefix", metavar="prefix"),
#   make_option(c("--numBoots"), type="numeric", default=NULL, 
#               help="Number of bootstraps to perform per individual", metavar="numeric")
# ); 

echo "Rscript $scriptDir/$script --infile $infile --outdir $indir --outPREFIX $outPREFIX --numBoots $numBoots"

Rscript $scriptDir/$script --infile $infile --outdir $outdir --outPREFIX $outPREFIX --numBoots $numBoots
done

sleep 10m
