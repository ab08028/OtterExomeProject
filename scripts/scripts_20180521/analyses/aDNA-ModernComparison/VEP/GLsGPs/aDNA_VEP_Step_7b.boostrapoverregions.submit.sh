#! /bin/bash
#$ -cwd
#$ -l h_rt=100:00:00,h_data=30G,highp
#$ -m abe
#$ -M ab08028
#$ -N bootstrapVEPResults
#$ -e /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -o /u/flashscratch/a/ab08028/captures/reports/angsd

source /u/local/Modules/default/init/modules.sh
module load R/3.5.1  



gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptDir=$gitDir/scripts/scripts_20180521/analyses/aDNA-ModernComparison/VEP/GLsGPs/
script=aDNA_VEP_Step_7b.boostrapoverregions.R

# parameters:
binsize=100000 # 100kb
numBoots=5 # eventually do more! 
minDepth=2 # I calculated totals with 1, 2 and 4. 2 lowers values of homAlt compared to 1. I think 2 is fitting since 1 read might seem homAlt but not be. 
minGP=0.95
minInd=1
dates='20190701-lowcov-AFprior-MajorMinor4 20190701-highcov-AFprior-MajorMinor4'
avgSitesFile=$SCRATCH/captures/aDNA-ModernComparison/VEP/sumGPsGLsPerVEPCategory/AVERAGECALLEDSITES.allInds.HighCov.LowCov.minDepth.${minDepth}.minInd.${minInd}.minGP.${minGP}.txt
# contains both dates
ref="mfur"
basename=angsdOut.mappedTo${ref}
type="GPs" # for now
# loop through dates:
for angsdDate in $dates
do
echo $angsdDate
indir=$SCRATCH/captures/aDNA-ModernComparison/angsd-GLs/$date/cdsPerCategoryFromVEP/
infile=$indir/${basename}.superfile.${type}.mafs.counts.0based.allCDSSites.AnnotatedWithVEP.bed.gz # need to fix header issue
chrSizes=$scriptDir/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.224.ChrLengths.txt # need the chr sizes to be uploaded somewhere
outdir=$SCRATCH/captures/aDNA-ModernComparison/compareMisSynDists_withBootstraps
mkdir -p $outdir
outPREFIX=${basename}.superfile.${type}
# file with avg sites to draw (depends on minDepth,minGP, minInd) and is a result of the script: /Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/analyses/aDNA-ModernComparison/VEP/GLsGPs/getAverageCalledCDSSItesAcrossAllInds.R
# get count for this date:
avgSitesToDraw=`grep $angsdDate $avgSitesFile | awk '{print $3}'` ## need to calculate! ##  can be found in /u/flashscratch/a/ab08028/captures/aDNA-ModernComparison/VEP/sumGPsGLsPerVEPCategory -- calculate average from across individuals (ask Kirk too)
echo "avg sites to draw: " $avgSitesToDraw

# parameters for R script:
# option_list = list(
#   make_option(c("-infile", "--infile"), type="character", default=NULL, 
#               help="path to your input file (superfile of cds sites with annotations from vep and GPs from angsd as well as depths", metavar="character"),
#   make_option(c("-chrSizes", "--chrSizes"), type="character", default=NULL, 
#               help="path to file with mustela chromosome sizes", metavar="character"),
#   make_option(c("-outdir", "--outdir"), type="character", default=NULL, 
#               help="path to outdir", metavar="character"),
#   make_option(c("-outPREFIX", "--outPREFIX"), type="character", default=NULL, 
#               help="outfilePrefix", metavar="character"),
#   make_option(c("-minDepth", "--minDepth"), type="character", default=NULL, 
#               help="min depth per individual for a site to be 'callalble'", metavar="character"),
#   make_option(c("-minGP", "--minGP"), type="character", default=NULL, 
#               help="minimum value of the max. GP per site, per individual. Use 0.95.", metavar="character"),
#   make_option(c("-binsize", "--binsize"), type="character", default=NULL, 
#               help="Size of bin to chunk the genome into (should be > than a recombination block)", metavar="character"),
#   make_option(c("-numBoots", "--numBoots"), type="character", default=NULL, 
#               help="Number of bootstraps to perform per individual", metavar="character"),
#   make_option(c("-avgSitesToDraw", "--avgSitesToDraw"), type="character", default=NULL, 
#               help="Number of sites to draw per individual (will be approximately this many, not exactly) due to varying bin size. Based on average 'callable' cds sites calculated elsewhere.", metavar="character")
# ); 
Rscript $scriptDir/$script --infile $infile --chrSizes $chrSizes --outdir $outdir --outPREFIX $outPREFIX --minDepth $minDepth --minGP $minGP --binsize $binsize --numBoots $numBoots --avgSitesToDraw $avgSitesToDraw
done