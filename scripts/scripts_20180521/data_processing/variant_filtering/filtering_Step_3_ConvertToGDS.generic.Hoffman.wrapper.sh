############## wrapper for GDS conversion script #############


module load R
gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/

scriptdir=${gitdir}/scripts/scripts_20180521/data_processing/variant_filtering/
script=filtering_Step_3_ConvertToGDS.generic.Hoffman.R

genotypeDate=20180806
vcfdir=/u/flashscratch/a/ab08028/captures/vcf_filtering/${genotypeDate}_filtered
outdir=$vcfdir/gdsFormat/
mkdir -p $outdir
# list of gzipped vcfs you want to convert (can be any VCFs)
vcfs2gds='snp_8_rmRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf.gz snp_8_rmRelativesAdmixed_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf.gz'
for vcf in $vcfs2gds
do
Rscript $scriptdir/${script} --vcf $vcfdir/$vcf --outdir $outdir
done
# output file name will be input file name with .gds extension instead of .vcf.gz extension
