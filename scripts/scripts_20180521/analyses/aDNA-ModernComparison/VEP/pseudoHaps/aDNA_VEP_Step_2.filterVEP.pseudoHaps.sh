#!/bin/bash
#$ -l h_rt=24:00:00,h_data=3G
#$ -N filter_vep
#$ -cwd
#$ -m bea
#$ -o /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -e /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -M ab08028

# Filter VEP output to separate missense, synonymous and LOF mutations that are in CANONICAL transcripts
# Then want to turn into coordinates 


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
#dates="20190524-highcov-AFprior 20190524-lowcov-AFprior 20190524-highcov-UNIFprior 20190524-lowcov-UNIFprior" # set of angsdDates you want to process 
dates="20190612-highcov-pseudoHaps 20190612-lowcov-pseudoHaps"
basename=angsdOut.mappedTo${ref}
for angsdDate in $dates
do
indir=$wd/VEP/pseudoHaps/$angsdDate 
input=${basename}.pseudoHaps.superfile.cdsOnly.1based.VEPInput.VEP.output.pick.tbl.gz

# 20190611: adding --pick to the initial run of VEP so that it only picks one annotation per variant
# this is what I thought filter-vep was doing, but on rare occasions it was still outputting a couple annotations per variant; don't want that. So for now, I'm just picking based on the VEP criteria (severity and validtity hierarchy). This doesn't change anything about this filtering script except the name of the input file. 

############## Filter Synonymous/NS ##########################
###### NS (missense):
$vepdir/filter_vep --filter "Consequence is missense_variant and CANONICAL is YES" \
--input_file $indir/$input \
--output_file $indir/${input%.tbl.gz}.missense.tbl \
--force_overwrite

###### LOF (stop_gained):
$vepdir/filter_vep --filter "Consequence is stop_gained and CANONICAL is YES" \
--input_file $indir/$input \
--output_file $indir/${input%.tbl.gz}.stopgained.tbl \
--force_overwrite

###### S (synonymous):
$vepdir/filter_vep --filter "Consequence is synonymous_variant and CANONICAL is YES" \
--input_file $indir/$input \
--output_file $indir/${input%.tbl.gz}.synonymous.tbl \
--force_overwrite

done

# can consider gzipping them 