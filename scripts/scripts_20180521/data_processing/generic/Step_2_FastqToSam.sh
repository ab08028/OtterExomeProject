#! /bin/bash
#$ -cwd
#$ -l h_rt=12:00:00,h_data=21G,arch=intel*
#$ -m bea

######## commenting out some of these settings so can reset when I submit the job: ###########

#### $ -o /u/flashscratch/a/ab08028/captures/reports
#### $ -e /u/flashscratch/a/ab08028/captures/reports
#### $ -M ab08028
########################################### 
# Usage: qsub -e $errorLocation -o $errorLocation -M $user -N fq2sam${c} $SCRIPTDIR/$scriptname [fastq R1] [fastqR2] [output bam name] [RGID : header_1a] [RGSM: sample ID] [RGLB: Lib1] [platform unit] [RGPL: platform] [seqCenter]

## This script will convert fastq file to uBAM file. You submit it to Hoffman using submit_run_step2


# RGID is unique identifier (you come up with, e.g. gidget_1a for gidget, library 1, lane a)
# RGSM is sample name (e.g. GIDGET)
# RGLB is library name (e.g. Lib1)
# PLATFORM_UNIT is a unique identifier for the sequencer found in the fastq file 
# RGPL is platform: illumina
# RGCN is sequencing center: UCB
# See the wrapper script, submit_run_FastqToSam.sh, for job submission example.

source /u/local/Modules/default/init/modules.sh
module load java/1.8.0_111 # need to be java 1.8

#PICARD=/u/home/a/ab08028/klohmueldata/annabel_data/bin/picard.jar
# update picard:
PICARD=/u/home/a/ab08028/klohmueldata/annabel_data/bin/Picard_2.8.1/picard.jar
SCRATCH=/u/flashscratch/a/ab08028
READ_DIR=$SCRATCH/captures/fastqs
BAM_OUTDIR=$SCRATCH/captures/bams
TEMP_DIR=$SCRATCH/temp
# note that read group ID should be unique for each sequenced sample; so same individual could
# be sequenced many times and have different read group IDs each time (e.g. CA_145_1a ; CA_145_1b)

mkdir -p ${BAM_OUTDIR}
java -Xmx16G -jar -Djava.io.tmpdir=$TEMP_DIR \
$PICARD FastqToSam \
FASTQ=$READ_DIR/$1 \
FASTQ2=$READ_DIR/$2 \
OUTPUT=$BAM_OUTDIR/$3 \
READ_GROUP_NAME=$4 \
SAMPLE_NAME=$5 \
LIBRARY_NAME=$6 \
PLATFORM_UNIT=$7 \
PLATFORM=$8 \
SEQUENCING_CENTER=$9 \
TMP_DIR=$TEMP_DIR

sleep 10m
