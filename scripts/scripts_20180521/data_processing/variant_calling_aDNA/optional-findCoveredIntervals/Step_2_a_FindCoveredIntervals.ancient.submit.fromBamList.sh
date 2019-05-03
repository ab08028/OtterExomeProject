#! /bin/bash
#$ -cwd
#$ -l h_rt=50:00:00,h_data=1G,highp
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
scriptDir=$gitDir/scripts/scripts_20180521/data_processing/variant_calling_aDNA
# script to run: 

#scriptname=Step_2_a_FindCoveredIntervals.ancient.sh # change this to final script name!! 
scriptname=Step_2_a_FindCoveredIntervals.ancient.fromBamList.sh 
# file locations:
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures
# for now just do ancient since they're done
#headers=$wd/samples/ancient.modern.comparison.txt # 3 ancient samples I'm processing further 
elutBamList=$scriptDir/angsd.bamList.mappedtoElutfullpaths.txt # list of bam files mapped to sea otter, including downsampled AND non-downsampled
mfurBamList=$scriptDir/angsd.bamList.mappedtoMfurfullpaths.txt  # list of bam files mapped to ferret, including downsampled AND non-downsampled

# references:
elutRef=/u/home/a/ab08028/klohmueldata/annabel_data/sea_otter_genome/dedup_99_indexed_USETHIS/sea_otter_23May2016_bS9RH.deduped.99.fasta
mfurRef=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta

reports=/u/flashscratch/a/ab08028/captures/reports/GATK/
mkdir -p $reports

# do elut first:

cat $elutBamList | while read bamPath
do
errorLocation=${reports}/angsd # report location
#mkdir -p $errorLocation
$QSUB -e $errorLocation -o $errorLocation -M $user -N intervals $scriptDir/$scriptname $bamPath $elutRef

#sleep 30s
done

# then do mfur
cat $mfurBamList | while read bamPath
do
errorLocation=${reports}/angsd # report location
#mkdir -p $errorLocation
$QSUB -e $errorLocation -o $errorLocation -M $user -N intervals $scriptDir/$scriptname $bamPath $mfurRef

#sleep 30s
done
