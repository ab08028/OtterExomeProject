#! /bin/bash
#$ -cwd
#$ -l h_rt=5:00:00,h_data=20G,arch=intel*
#$ -m bea

QSUB=/u/systems/UGE8.0.1vm/bin/lx-amd64/qsub

# location of github:  which may be on remote server or local laptop
gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
# scripts:
scriptDir=$gitDir/scripts/scripts_20180521/data_processing/
# script to run: 
scriptname=Step_1_b_AdapterRemoval.sh

# file locations:
SCRATCH=/u/flashscratch/a/ab08028/
wd=$SCRATCH/captures
fastqs=$wd/fastqs
adapterRemoval=$wd/adapterRemoval
mkdir -p $wd/adapterRemoval
bams=$wd/bams

# job info: 
errorLocation=/u/flashscratch/a/ab08028/captures/reports/step_1_b_adapterRemoval # report location
user=ab08028 # where emails are sent

## usage: qsub -e $errorLocation -o $errorLocation -M $user -N adaptRemoval${c} $scriptDir/generic/$scriptname $fastqs/$fileR1 $fastqs/$fileR2
### Adapter removal with AdapterRemoval2.0
# paired end mode; collapse reads


# starting prefix number
START=1
# ending prefix number 
END=12
# ancient DNA is treated separately because start with A1... # 
cd $fastqs
for (( c=$START; c<=$END; c++ ))
do
fileR1=`ls A${c}_Elut*R1*fastq.gz` # the R1 fastq file
fileR2=`ls A${c}_Elut*R2*fastq.gz` # the R2 fastq file
header=${fileR1%_S*_R*} # this is the header sample name
# usage: qsub -e $errorLocation -o $errorLocation -M $user -N adaptRemov${c} $scriptDir/generic/$scriptname $fastqs/[fastq input R1] $fastqs/[fastq input R2]  [basename header for output]
$QSUB -e $errorLocation -o $errorLocation -M $user -N adaptRemov${c} $scriptDir/generic/$scriptname $fastqs/$fileR1 $fastqs/$fileR2 $wd/adapterRemoval/${header}
# clear variables just in case:
fileR1=""
fileR2=""
header=""
sleep 10m # sleep between job submissions
done
