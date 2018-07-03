#! /bin/bash
#$ -cwd
#$ -l h_rt=300:00:00,h_data=25G,arch=intel*,highp
#$ -t 8
#$ -m bea

####### run paleomix

#### load modules: ####
source /u/local/Modules/default/init/modules.sh
module load python/2.7
module load java/1.8.0_111
module load samtools
module load bwa
module load R
# need mapdamage in the path! 
######## directories #######
SCRATCH=/u/flashscratch/a/ab08028

######## run paleomix: ############
###### BEFORE YOU RUN PALEOMIX, do a dry run: (replace run with dry_run) #### 
# this will make sure everything you need is in place and will make indices for you
# Do this once before you set the scripts in motion
plmx=/u/home/a/ab08028/klohmueldata/annabel_data/bin/paleomixLinks/paleomix
MAX_THREADS=8 # adjust as needed 
makefile=$1
DESTINATION=$2 # location where you want files to go 
TEMP_DIR=$SCRATCH/temp
JAR_ROOT=/u/home/a/ab08028/klohmueldata/annabel_data/bin/Picard_2.8.1/
$plmx bam_pipeline run $makefile \
--max-threads=${MAX_THREADS} \
--adapterremoval-max-threads=${MAX_THREADS} \
--bwa-max-threads=${MAX_THREADS} \
--temp-root=${TEMP_DIR} \
--destination=${DESTINATION} \
--jar-root=${JAR_ROOT}

# test: 
