#!/bin/bash
#$ -l h_rt=40:00:00,h_data=3G,highp
#$ -N aDNA_VEP
#$ -pe shared 3
#$ -cwd
#$ -m bea
#$ -o /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -e /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -M ab08028
################# Run VEP on cds input from angsd #############
### then have to figure out a way to 
# to install:
# module load perl
# module load htslib/1.3.2
# run with --NO_HTSLIB
# say y to cache file; install 88 (mustela)
######### RUNNING ON CDS SITES ONLY ##########
# load your modules:
source /u/local/Modules/default/init/modules.sh
module load perl/5.10.1
module load htslib
vepdir=/u/home/a/ab08028/klohmueldata/annabel_data/bin/ensembl-vep/

# directories: 
gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptDir=$gitDir/scripts/scripts_20180521/

SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/aDNA-ModernComparison
#GLdir=$wd/angsd-GLs

ref=mfur # only work with MFUR when going into VEP!
dates="20190524-highcov-AFprior 20190524-lowcov-AFprior 20190524-highcov-UNIFprior 20190524-lowcov-UNIFprior" # set of angsdDates you want to process 
basename=angsdOut.mappedTo${ref}
for angsdDate in $dates
do
indir=$wd/VEP/$angsdDate 

# can use GP or GL superfile; they are in the same order and have same sites, so it doesn't matter, since GLs and GPs aren't taken along for the ride in the VEP input format
# I am using GPs, but again, could be either:
vepinput=${basename}.superfile.GPs.mafs.counts.cdsOnly.1based.VEPInput.txt.gz # may not work when gzipped -- expt and see.

### adding CANONICAL field so I can filter on that 
$vepdir/vep -v -i ${indir}/$vepinput --fork 3 \
--cache --force_overwrite --species mustela_putorius_furo \
--numbers --domains --variant_class --canonical \
-o $indir/${vepinput%.txt.gz}.VEP.output.tbl

# gzip output:
gzip -f $indir/${vepinput%.txt.gz}.VEP.output.tbl

# note if you use a gzipped file you'll get the error "gzip: stdout: Broken pipe" but it doesn't actually break anything
done
### will have to figure out how to merge this back with my angsd cds super file (maybe based on the "marker" column>)
# MAKE SURE TO DEAL WITH NON UNIQUE ENTRIES CAREFULLY
# LOOK AT MY OLD FILTERING VEP RESULTS SCRIPTS FOR INSPIRATION : /Users/annabelbeichman/Documents/UCLA/Otters/otterScriptsGithub/OtterGenomeProject/DiversityAnalyses/northernSeaOtterAnalyses/VEP/VEP.GeneticLoad.Scripts
#sleep 5m
