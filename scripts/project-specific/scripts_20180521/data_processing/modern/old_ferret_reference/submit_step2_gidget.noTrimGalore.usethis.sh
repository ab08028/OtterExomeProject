#! /bin/bash
#$ -cwd
#$ -l h_rt=24:00:00,h_data=500M
#$ -N subFq2sam
#$ -o /u/scratch/a/ab08028/otters/reports/
#$ -e /u/scratch/a/ab08028/otters/reports/
#$ -m bea
#$ -M ab08028

# this script will submit run_step2_FastqToSam.AB.sh that converts fastq to uBAM file format

SCRIPTDIR=/u/home/a/ab08028/klohmueldata/annabel_data/mapReads_General/scripts

QSUB=/u/systems/UGE8.0.1vm/bin/lx-amd64/qsub

# December 2016: mapping 2 lanes of Gidget shotgun sequencing (HiSeq4000 150 PE) 
# job name , script, input fq 1, input fq 2, output bam name, read group ID (sample_libLane), sample name, library name, flow cell identifier, sequencing center (UCB = Berkeley)
## 20170630: Just noticed: I didn't map the "val" fastqs for Gidget, I did the unfiltered fastqs. Need to redo. So add Gidget here.

## GIDGET:
# command I used : this was actually CORRECT! According to JAR's NEW pipeline, you don't do TrimGalore -- use MarkIlluminaAdapters instead!
$QSUB -N fq2sam01 $SCRIPTDIR/run_step2_FastqToSam.AB.sh RWAB001_L002_R1_001.fastq.gz RWAB001_L002_R2_001.fastq.gz RWAB001_L002_FastqToSam.bam GIDGET_1a GIDGET Lib1 H5FLGBBXX.2 UCB
# why was it  Step_2_FastqToSam.sh
run_step2_FastqToSam.AB.sh
# don't do this on post trim-galore reads: $QSUB -N fq2sam01 $SCRIPTDIR/run_step2_FastqToSam.AB.sh RWAB001_L002_R1_001_val_1.fq.gz RWAB001_L002_R2_001_val_2.fq.gz RWAB001_L002_FastqToSam.bam GIDGET_1a GIDGET Lib1 H5FLGBBXX.2 UCB

sleep 30m
# command I used : this was actually CORRECT! According to JAR's NEW pipeline, you don't do TrimGalore -- use MarkIlluminaAdapters instead!
$QSUB -N fq2sam02 $SCRIPTDIR/run_step2_FastqToSam.AB.sh RWAB001_S6_L006_R1_001.fastq.gz  RWAB001_S6_L006_R2_001.fastq.gz RWAB001_S6_L006_FastqToSam.bam GIDGET_1b GIDGET Lib1 H7MWWBBXX.6 UCB


# don't do this on post trim-galore reads: $QSUB -N fq2sam01 $SCRIPTDIR/run_step2_FastqToSam.AB.sh RWAB001_S6_L006_R1_001_val_1.fq.gz RWAB001_S6_L006_R2_001_val_2.fq.gz RWAB001_L002_FastqToSam.bam GIDGET_1b GIDGET Lib1 H5FLGBBXX.2 UCB


