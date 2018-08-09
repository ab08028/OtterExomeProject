#! /bin/bash
#$ -cwd
#$ -l h_rt=00:10:00,h_data=16G
#$ -N countNoCall_perInd
#$ -o /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -e /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -m abe
#$ -M ab08028
source /u/local/Modules/default/init/modules.sh
module load python/2.7
#### file locations
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/vcf_filtering
scriptdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_filtering
script=filtering_getNoCallPerInd.py
rundate=20180724 # date genotypes were called
vcfdir=$wd/${rundate}_filtered
outdir=$vcfdir/filteringStats
mkdir -p $outdir
python $scriptdir/$script $vcfdir/dummyVCF.forSandbox.allSites_5_passingFilters.vcf.gz $outdir/dummy.filtering.test.PerInd.txt

