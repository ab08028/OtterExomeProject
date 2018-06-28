#! /bin/bash
#$ -wd /u/scratch/a/ab08028/otters/bams
#$ -l h_rt=24:00:00,h_data=1G
#$ -N subBamFltr
#$ -o /u/scratch/a/ab08028/otters/reports
#$ -e /u/scratch/a/ab08028/otters/reports
#$ -m abe
#$ -M ab08028
# do this after indel realignment (reduces bam size a lot)
# do qualimap before and after
SCRIPT_DIR=/u/home/a/ab08028/klohmueldata/annabel_data/mapReads_General/scripts
QSUB=/u/systems/UGE8.0.1vm/bin/lx-amd64/qsub
INDIR=/u/scratch2/a/ab08028/otters/bams/IndelRealigned

#cd /u/scratch/j/jarobins/irnp/bams/MarkDup/done
# cd /u/scratch/a/ab08028/otters/bams/MarkDup
# for i in {01..19}; do
# 	FILE=${i}*bam
# 	${QSUB} -N bamfiltr${i} $SCRIPT_DIR/run_step6_RemoveBadReads_110716.sh $INDIR/${FILE}
# 	sleep 1h
# done


# manually for Gidget:
qsub /u/home/a/ab08028/klohmueldata/annabel_data/mapReads_General/scripts/run_step7_filterBadReads-altNoIndelRealignment.optional.sh 01_Elut_CA_Gidget_Aligned_MarkDup.bam