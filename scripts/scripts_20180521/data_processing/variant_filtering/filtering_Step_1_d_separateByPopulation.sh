#! /bin/bash
#$ -cwd
#$ -l h_rt=50:00:00,h_data=16G,highp,arch=intel*
#$ -N 1b_vcf_filtering
#$ -o /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -e /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -m abe
#$ -M ab08028

rundate=20180806
#### file locations
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/vcf_filtering
infile=raw_variants.vcf.gz ### make sure this doesn't have a path as part of its name! just infile names
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta
badIndFile=$wd/bad.individuals.toRemove.txt # IDs of individuals to remove (get from step 1a)
# need to remove bad individuals from lists somehow 
popHeaderDir=$SCRATCH/captures/samples/popHeaders_rmBadInd/


#################################################################################
############################ Separate by population #############################
#################################################################################
# want to separate vcfs by population