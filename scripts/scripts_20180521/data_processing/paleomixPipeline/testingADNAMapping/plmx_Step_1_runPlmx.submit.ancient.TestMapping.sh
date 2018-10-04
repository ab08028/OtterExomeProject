#! /bin/bash
#$ -cwd
#$ -l h_rt=50:00:00,h_data=1G,highp
#$ -o /u/flashscratch/a/ab08028/captures/reports/submissions/
#$ -e /u/flashscratch/a/ab08028/captures/reports/submissions/
#$ -m bea
#$ -M ab08028
#$ -N plmx_submit
#***************** REMEMBER TO SUBMIT with qsub -V !!!! ***********
######### This script run will submit a series of jobs that convert fastq to sam and adds readgroup info
QSUB=/u/systems/UGE8.0.1vm/bin/lx-amd64/qsub
# location of github:  which may be on remote server or local laptop
gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
# scripts:
scriptDir=$gitDir/scripts/scripts_20180521/data_processing/paleomixPipeline
# script to run: 
scriptname=plmx_Step_1_runPlmx.sh

# file locations:
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures
makefileDir=$scriptDir/makefiles/ancientMakefiles-TestMapping
#headers=$wd/samples/ancientSamples.txt
#headers=$wd/samples/ancientSamples.temp.remainder.txt
headers=$wd/samples/ancientLibs.Screen.txt # 20181004: screening libraries
# outdirectory:
outdir=$wd/paleomix/testMapping
mkdir -p $outdir
# job info: 
user=ab08028 # where emails are sent

# usage; qsub script [makefile, full path] [outdir]

cat $headers | while read header
do
errorLocation=/u/flashscratch/a/ab08028/captures/reports/paleomix/testMapping/${header} # report location
mkdir -p $errorLocation
mkdir -p $outdir/${header}
$QSUB -e $errorLocation -o $errorLocation -M $user -N plmx.${header} \
$scriptDir/$scriptname $makefileDir/${header}.paleomix.makefile.yaml $outdir/${header}
# clear variables:
errorLocation=""
sleep 5m
done


