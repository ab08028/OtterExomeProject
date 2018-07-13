#! /bin/bash
#$ -cwd
<<<<<<< HEAD
#$ -l h_rt=50:00:00,h_data=1G,highp
=======
#$ -l h_rt=5:00:00,h_data=1G,highp
>>>>>>> aa608cd3e5b551aac391c322e90cf4610353b5ed
#$ -o /u/flashscratch/a/ab08028/captures/reports/submissions/
#$ -e /u/flashscratch/a/ab08028/captures/reports/submissions/
#$ -m bea
#$ -M ab08028
#$ -N intervals_submit
user=ab08028 # where emails are sent

QSUB=/u/systems/UGE8.0.1vm/bin/lx-amd64/qsub

# location of github:  which may be on remote server or local laptop
gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
# scripts:
scriptDir=$gitDir/scripts/scripts_20180521/data_processing/variant_calling
# script to run: 
<<<<<<< HEAD
scriptname=Step_2_a_FindCoveredIntervals.sh # change this to final script name!! 
=======
scriptname=Step_2_a_FindCoveredIntervals.ancient.sh # change this to final script name!! 
>>>>>>> aa608cd3e5b551aac391c322e90cf4610353b5ed
# 
# file locations:
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures
<<<<<<< HEAD
# for now just do ancient since they're done
headers=$wd/samples/ancientSamples.txt # all ancient (for now; waiting on modern to finish)
=======
headers=$wd/samples/ancientSamples.txt # all Elut samples, modern, ancient and blank
>>>>>>> aa608cd3e5b551aac391c322e90cf4610353b5ed
reports=/u/flashscratch/a/ab08028/captures/reports/GATK/
mkdir -p $reports

cat $headers | while read header
do
errorLocation=${reports}/${header} # report location
mkdir -p $errorLocation
$QSUB -e $errorLocation -o $errorLocation -M $user -N intervals.${header} $scriptDir/$scriptname $header
<<<<<<< HEAD
sleep 10m
=======
sleep 30s
>>>>>>> aa608cd3e5b551aac391c322e90cf4610353b5ed
done
