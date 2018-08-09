######### wrapper for the get no call count per individual python script
#### file locations
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/vcf_filtering
scriptdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_filtering
script=filtering_getNoCallPerInd.py
#### parameters:
rundate=20180724 # date genotypes were called

vcfdir=$wd/${rundate}_filtered
outdir=$vcfdir/filterStats
mkdir -p $outdir
# this could be any vcf you want. maybe all of them? 
python $scriptdir/$script $vcfdir/all_6_passingBespoke_passingFilters_80percCall_raw_variants.vcf.gz $outdir/all_6_passingBespoke_passingFilters_80percCall_raw_variants.NoCall.PerInd.txt