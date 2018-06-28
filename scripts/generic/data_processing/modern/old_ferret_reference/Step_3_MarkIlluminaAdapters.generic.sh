#! /bin/bash
#$ -wd /u/scratch/a/ab08028/otters/
#$ -l h_rt=12:00:00,h_data=21G,arch=intel*,highp
#$ -o /u/scratch/a/ab08028/otters/reports
#$ -e /u/scratch/a/ab08028/otters/reports
#$ -m bea
#$ -M ab08028

# Usage: qsub -N jobname run_step3_MarkIlluminaAdapters.AB.sh INPUT_FastqToSam.bam

# why is this an array?
source /u/local/Modules/default/init/modules.sh
module load java

wd=workingDirectory # set working dir

cd $wd/bams

FILENAME=$1
PICARD=picardLocation

TEMP_DIR=$wd/temp

java -Xmx16G -jar -Djava.io.tmpdir=$TEMP_DIR \
$PICARD \
MarkIlluminaAdapters \
I=$FILENAME \
O=${FILENAME%_FastqToSam.bam}_MarkIlluminaAdapters.bam \
M=${FILENAME%_FastqToSam.bam}_MarkIlluminaAdapters.bam_metrics.txt \
TMP_DIR=$TEMP_DIR
