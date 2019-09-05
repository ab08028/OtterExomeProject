#! /bin/bash
#$ -cwd
#$ -l h_rt=100:00:00,h_data=30G,highp
#$ -m abe
#$ -M ab08028
#$ -N setUpBins
#$ -e /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -o /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -V

source /u/local/Modules/default/init/modules.sh
module load R/3.5.1  

### each site should only appear once in the super file but want to make sure 
# zcat angsdOut.mappedTomfur.superfile.GPs.mafs.counts.0based.allCDSSites.AnnotatedWithVEP.bed.gz | grep -v "#" -c   # total sites
# 25971997
# zcat angsdOut.mappedTomfur.superfile.GPs.mafs.counts.0based.allCDSSites.AnnotatedWithVEP.bed.gz | grep -v "#" | awk '{print $4}' | sort | uniq | wc -l
# 25971997 with uniq
# so since --pick was used, each site only appears once
SCRATCH=/u/flashscratch/a/ab08028/

gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptDir=$gitDir/scripts/scripts_20180521/analyses/aDNA-ModernComparison/VEP/GLsGPs/
#script=NEW.aDNA_VEP_Step_7b-i-alt.setUpBootstrapBinsPerAncientModern.R
script=NEW.aDNA_VEP_Step_7b-i-alt-2.setUpBootstrapBinsPerPopulationAndAncientModern.R 
# parameters:
binsize=500000 # 500kb (raised from 100kb)
#numBoots=100 # eventually do more! 
minDepth=2 # I calculated totals with 1, 2 and 4. 2 lowers values of homAlt compared to 1. I think 2 is fitting since 1 read might seem homAlt but not be. 
minGP=0.95
minInd=1
minCalledSitesPerWindowPerIndividual=1000 # for now
minIndPerWindow=9 # all inds must be represented in the window, though not necessarily overlapping
type="GPs" 
ref="mfur"
################## low coverage #####################
angsdDate="20190701-lowcov-AFprior-MajorMinor4"
echo $angsdDate
bamList="/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_calling_aDNA/bamLists/SampleIDsInOrder.LowCoverageOnly.BeCarefulOfOrder.txt"
# contains both dates
basename=angsdOut.mappedTo${ref}


outPREFIX=${basename}.Bins.${type}.ProbCutoff.${minGP}.DepthCutoff.${minDepth}.minInd.${minInd}.${angsdDate}

indir=$SCRATCH/captures/aDNA-ModernComparison/angsd-GLs/$angsdDate/cdsPerCategoryFromVEP/
#infile=$indir/testingScript.CDS.bed.gz
infile=$indir/${basename}.superfile.${type}.mafs.counts.0based.allCDSSites.AnnotatedWithVEP.bed.gz # need to fix header issue
chrSizes=$scriptDir/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.224.ChrLengths.txt # need the chr sizes to be uploaded somewhere
outdir=$SCRATCH/captures/aDNA-ModernComparison/VEP/compareMisSynDists_withBootstraps/separateCA-AK/$angsdDate/
mkdir -p $outdir


# option_list = list(
#   make_option(c("--infile"), type="character", default=NULL, 
#               help="path to your input file (superfile of cds sites with annotations from vep and GPs from angsd as well as depths", metavar="file"),
#   make_option(c("--chrSizes"), type="character", default=NULL, 
#               help="path to file with mustela chromosome sizes", metavar="character"),
#   make_option(c("--outdir"), type="character", default=NULL, 
#               help="path to outdir", metavar="dir"),
#   make_option(c("--outPREFIX"), type="character", default=NULL, 
#               help="outfilePrefix", metavar="prefix"),
#   make_option(c("--minDepth"), type="numeric", default=NULL, 
#               help="min depth per individual for a site to be 'callalble'", metavar="numeric"),
#   make_option(c("--minGP"), type="numeric", default=NULL, 
#               help="minimum value of the max. GP per site, per individual. Use 0.95.", metavar="numeric"),
#   make_option(c("--minCalledSitesPerWindowPerIndividual"), type="numeric", default=NULL, 
#               help="minimum sites called per individual per window", metavar="numeric"),
#   make_option(c("--minIndPerWindow"), type="numeric", default=NULL, 
#               help="minimum individuals with data per window; data does not need to be overlapping, but most contain at least minCalledSitesPerWindowPerIndividual in the window", metavar="numeric"),
#   make_option(c("--binsize"), type="numeric", default=NULL, 
#               help="Size of bin to chunk the genome into (should be > than a recombination block)", metavar="numeric"),
#   make_option(c("--bamList"), type="character", default=NULL, 
#               help="path to list of bams in angsd (gives IDs) ***ASSUMES THAT ANCIENT SAMPLE IDs start with 'A'!!!!!!#", metavar="file")
# ); 

echo "Rscript $scriptDir/$script --infile $infile --chrSizes $chrSizes --outdir $outdir --outPREFIX $outPREFIX --minDepth $minDepth --minGP $minGP --binsize $binsize --minCalledSitesPerWindowPerIndividual $minCalledSitesPerWindowPerIndividual --minIndPerWindow $minIndPerWindow --bamList $bamList"

Rscript $scriptDir/$script --infile $infile --chrSizes $chrSizes --outdir $outdir --outPREFIX $outPREFIX --minDepth $minDepth --minGP $minGP --binsize $binsize --minCalledSitesPerWindowPerIndividual $minCalledSitesPerWindowPerIndividual --minIndPerWindow $minIndPerWindow --bamList $bamList

################## low coverage #####################
angsdDate="20190701-highcov-AFprior-MajorMinor4"
echo $angsdDate
bamList="/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_calling_aDNA/bamLists/SampleIDsInOrder.HighCoverageAndADNAOnly.BeCarefulOfOrder.txt"
# contains both dates
basename=angsdOut.mappedTo${ref}

# loop through dates:
outPREFIX=${basename}.Bins.${type}.ProbCutoff.${minGP}.DepthCutoff.${minDepth}.minInd.${minInd}.${angsdDate}

indir=$SCRATCH/captures/aDNA-ModernComparison/angsd-GLs/$angsdDate/cdsPerCategoryFromVEP/
#infile=$indir/testingScript.CDS.bed.gz
infile=$indir/${basename}.superfile.${type}.mafs.counts.0based.allCDSSites.AnnotatedWithVEP.bed.gz # need to fix header issue
chrSizes=$scriptDir/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.224.ChrLengths.txt # need the chr sizes to be uploaded somewhere
outdir=$SCRATCH/captures/aDNA-ModernComparison/VEP/compareMisSynDists_withBootstraps/separateCA-AK/$angsdDate/
mkdir -p $outdir

# option_list = list(
#   make_option(c("--infile"), type="character", default=NULL, 
#               help="path to your input file (superfile of cds sites with annotations from vep and GPs from angsd as well as depths", metavar="file"),
#   make_option(c("--chrSizes"), type="character", default=NULL, 
#               help="path to file with mustela chromosome sizes", metavar="character"),
#   make_option(c("--outdir"), type="character", default=NULL, 
#               help="path to outdir", metavar="dir"),
#   make_option(c("--outPREFIX"), type="character", default=NULL, 
#               help="outfilePrefix", metavar="prefix"),
#   make_option(c("--minDepth"), type="numeric", default=NULL, 
#               help="min depth per individual for a site to be 'callalble'", metavar="numeric"),
#   make_option(c("--minGP"), type="numeric", default=NULL, 
#               help="minimum value of the max. GP per site, per individual. Use 0.95.", metavar="numeric"),
#   make_option(c("--minCalledSitesPerWindowPerIndividual"), type="numeric", default=NULL, 
#               help="minimum sites called per individual per window", metavar="numeric"),
#   make_option(c("--minIndPerWindow"), type="numeric", default=NULL, 
#               help="minimum individuals with data per window; data does not need to be overlapping, but most contain at least minCalledSitesPerWindowPerIndividual in the window", metavar="numeric"),
#   make_option(c("--binsize"), type="numeric", default=NULL, 
#               help="Size of bin to chunk the genome into (should be > than a recombination block)", metavar="numeric"),
#   make_option(c("--bamList"), type="character", default=NULL, 
#               help="path to list of bams in angsd (gives IDs) ***ASSUMES THAT ANCIENT SAMPLE IDs start with 'A'!!!!!!#", metavar="file")
# ); 

echo "Rscript $scriptDir/$script --infile $infile --chrSizes $chrSizes --outdir $outdir --outPREFIX $outPREFIX --minDepth $minDepth --minGP $minGP --binsize $binsize --minCalledSitesPerWindowPerIndividual $minCalledSitesPerWindowPerIndividual --minIndPerWindow $minIndPerWindow --bamList $bamList"

Rscript $scriptDir/$script --infile $infile --chrSizes $chrSizes --outdir $outdir --outPREFIX $outPREFIX --minDepth $minDepth --minGP $minGP --binsize $binsize --minCalledSitesPerWindowPerIndividual $minCalledSitesPerWindowPerIndividual --minIndPerWindow $minIndPerWindow --bamList $bamList


sleep 5m
