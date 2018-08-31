#! /bin/bash
#$ -cwd
#$ -l h_rt=24:00:00,h_data=16G
<<<<<<< HEAD
#$ -N bespoke_vcf_filtering
=======
#$ -N countNoCall_perInd
>>>>>>> 8fa0029804ee0f1f6d26b1d31106eed5eb5dac33
#$ -o /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -e /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -m abe
#$ -M ab08028
<<<<<<< HEAD
=======

>>>>>>> 8fa0029804ee0f1f6d26b1d31106eed5eb5dac33
######### wrapper for the get no call count per individual python script
## need to use python 2.7 or later!
source /u/local/Modules/default/init/modules.sh
module load python/2.7
#### file locations
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/vcf_filtering
scriptdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_filtering
script=filtering_getNoCallPerInd.py
#### parameters:
rundate=20180724 # date genotypes were called

vcfdir=$wd/${rundate}_filtered
outdir=$vcfdir/filteringStats
mkdir -p $outdir
# this could be any vcf you want. maybe all of them? 
python $scriptdir/$script $vcfdir/all_6_passingBespoke_passingFilters_80percCall_raw_variants.vcf.gz $outdir/all_6_passingBespoke_passingFilters_80percCall_raw_variants.NoCall.PerInd.txt


# this could be any vcf you want. maybe all of them? 
<<<<<<< HEAD
python $scriptdir/$script $vcfdir/all_6_passingBespoke_passingFilters_80percCall_raw_variants.vcf.gz $outdir/all_6_passingBespoke_passingFilters_80percCall_raw_variants.NoCall.PerInd.txt
=======
#python $scriptdir/$script $vcfdir/all_6_passingBespoke_passingFilters_80percCall_raw_variants.vcf.gz $outdir/all_6_passingBespoke_passingFilters_80percCall_raw_variants.NoCall.PerInd.txt
>>>>>>> 4b497b1d52d2c08ed41c2cdfb49756adc09bc2e4
