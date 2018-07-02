#! /bin/bash
#$ -cwd
#$ -l h_rt=300:00:00,h_data=25G,arch=intel*
#$ -t 8
#$ -m bea

####### run paleomix

#### load modules: ####
source /u/local/Modules/default/init/modules.sh
module load python/2.7
module load java/1.8.0_111

######## directories #######
SCRATCH=/u/flashscratch/a/ab08028

######## run paleomix: ############

plmx=/u/home/a/ab08028/klohmueldata/annabel_data/bin/paleomixLinks/paleomix
MAX_THREADS=8 # adjust as needed 
makefile=$1
DESTINATION=$2 # location where you want files to go 
TEMP_DIR=$SCRATCH/temp
$plmx bam_pipeline run $makefile \
--max-threads=${MAX_THREADS} \
--adapterremoval-max-threads=${MAX_THREADS} \
--bwa-max-threads=${MAX_THREADS} \
--temp-root=${TEMP_DIR} \
--destination=${DESTINATION}

# test: 
