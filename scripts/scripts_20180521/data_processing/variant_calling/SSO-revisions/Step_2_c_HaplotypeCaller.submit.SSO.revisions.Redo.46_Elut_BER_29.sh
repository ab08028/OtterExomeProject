#! /bin/bash
#$ -cwd
#$ -l h_rt=3:00:00,h_data=3G,highp
#$ -o /u/flashscratch/a/ab08028/captures/reports/submissions/
#$ -e /u/flashscratch/a/ab08028/captures/reports/submissions/
#$ -m bea
#$ -M ab08028
#$ -N hapCaller_submit_SSO
user=ab08028 # where emails are sent

QSUB=/u/systems/UGE8.6.4/bin/lx-amd64/qsub # trying updated version here 
REFSHORTCODE=SSO 
# location of github:  which may be on remote server or local laptop
gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
# scripts:
scriptDir=$gitDir/scripts/scripts_20180521/data_processing/variant_calling/${REFSHORTCODE}-revisions/

# script to run: 
scriptname=Step_2_c_HaplotypeCaller.${REFSHORTCODE}.revisions.sh # change this to final script name!! 
# 
# file locations:
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures
#headers=$wd/samples/SamplesPartOf.20181119.JointGenotyping.UsedForMappingInRevision.includesRWAB.txt # note headers
header=46_Elut_BER_29
reports=/u/flashscratch/a/ab08028/captures/reports/GATK/
mkdir -p $reports

#cat $headers | while read header
#do
errorLocation=${reports}/${header} # report location
mkdir -p $errorLocation
$QSUB -e $errorLocation -o $errorLocation -M $user -N hapCaller.${REFSHORTCODE}.${header} $scriptDir/$scriptname $header
sleep 1m
#done
