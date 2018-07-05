#! /bin/bash
#$ -cwd
#$ -l h_rt=5:00:00,h_data=1G,highp
#$ -o /u/flashscratch/a/ab08028/captures/reports/submissions/
#$ -e /u/flashscratch/a/ab08028/captures/reports/submissions/
#$ -m bea
#$ -M ab08028
#$ -N fastqc_submit
user=ab08028 # where emails are sent

QSUB=/u/systems/UGE8.0.1vm/bin/lx-amd64/qsub

# location of github:  which may be on remote server or local laptop
gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
# scripts:
scriptDir=$gitDir/scripts/scripts_20180521/data_processing/variant_calling
# script to run: 
scriptname=Step_0_c_FastQC.Hoffman.sh # change this to final script name!! 
# 
# file locations:
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures
headers=$wd/samples/allElutSamples.txt # all Elut samples, modern, ancient and blank
errorLocation=/u/flashscratch/a/ab08028/captures/reports/fastqc/
mkdir -p $errorLocation

cat $headers | while read header
do
$QSUB -e $errorLocation -o $errorLocation -M $user -N fastqc.${header} $scriptDir/$scriptname $header
sleep 10s
done