#! /bin/bash
#$ -wd /u/scratch/a/ab08028/otters/bams/MarkDup
#$ -l h_rt=14:00:00,h_data=1G
#$ -N subMrkDup
#$ -o /u/scratch/a/ab08028/otters/bams/MarkDup
#$ -e /u/scratch/a/ab08028/otters/bams/MarkDup
#$ -m abe
#$ -M ab08028

SCRIPT_DIR=/u/home/a/ab08028/klohmueldata/annabel_data/mapReads_General/scripts
QSUB=/u/systems/UGE8.0.1vm/bin/lx-amd64/qsub

$QSUB -N MarkDup04 $SCRIPTDIR/run_step5_MarkDuplicates.AB.sh 04_rwjr002_S2_L002_Aligned.bam 03_IRNP_CL065_Aligned_MarkDup.bam
sleep 2h
$QSUB -N MarkDup05 $SCRIPTDIR/run_MarkDuplicates.sh 05_rwjr003_S3_L003_Aligned.bam 04_IRNP_RKW2524_Aligned_MarkDup.bam
sleep 2h
$QSUB -N MarkDup06 $SCRIPTDIR/run_MarkDuplicates.sh 06_rwjr004_S4_L004_Aligned.bam 05_IRNP_CL075_Aligned_MarkDup.bam
sleep 2h
$QSUB -N MarkDup07 $SCRIPTDIR/run_MarkDuplicates.sh 07_rwjr005_S5_L005_Aligned.bam 06_IRNP_CL061_Aligned_MarkDup.bam
sleep 2h
$QSUB -N MarkDup08 $SCRIPTDIR/run_MarkDuplicates.sh 08_rwjr006_S6_L006_Aligned.bam 07_IRNP_CL149_Aligned_MarkDup.bam
sleep 2h
$QSUB -N MarkDup09 $SCRIPTDIR/run_MarkDuplicates.sh 09_rwjr007_S7_L007_Aligned.bam 08_IRNP_RKW119_Aligned_MarkDup.bam
sleep 2h
$QSUB -N MarkDup10 $SCRIPTDIR/run_MarkDuplicates.sh 10_rwjr008_S8_L008_Aligned.bam 09_IRNP_CL189_Aligned_MarkDup.bam
