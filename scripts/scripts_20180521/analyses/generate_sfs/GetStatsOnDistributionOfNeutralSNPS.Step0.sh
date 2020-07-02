####### first want to exclude monomorphic 
module load vcftools
wd=$SCRATCH/captures/basicSNPstats/distOfNeutralSNPs
mkdir -p $wd
cd $wd
vcfdir=/u/flashscratch/a/ab08028/captures/vcf_filtering/20181119_filtered/neutral_and_cds_VCFs/neutralVCFs
infile=neutral.snp_9b_forEasySFS_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf
pop=AK
#outfile=$wd/${infile%.vcf.gz}.mac1Filter.NoMonomorphic.ForGettingDistOfSNPs.$pop
# want to do mac 1 filter and filter for high degrees of missing data
#vcftools --gzvcf $infile --mac 1 --recode --out $outfile 
# --max-missing 0.2


outfile=$pop.intFile1
# get pop specific vcf (with mac1):
vcftools --gzvcf $vcfdir/$infile --recode --out $wd/$outfile  --mac 1 --keep /u/flashscratch/a/ab08028/captures/samples/popsForFst/$pop.forFst.txt
# then want to filter within population:
vcftools --vcf $wd/$pop.intFile1.recode.vcf --recode --mac 1 --max-missing-count 4 --out $wd/$pop.Mac1.maxMiss4