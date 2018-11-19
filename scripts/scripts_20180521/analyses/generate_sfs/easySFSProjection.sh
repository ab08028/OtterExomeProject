####### Easy SFS
# https://github.com/isaacovercast/easySFS
# install:
# git clone git@github.com:isaacovercast/easySFS.git
# cd easySFS
# chmod +x *.py
# easySFS.py

vcfdir=/u/flashscratch/a/ab08028/captures/vcf_filtering/20180806_filtered/
vcf=all_8_rmRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_0.9_rmBadIndividuals_passingFilters_raw_variants.vcf.gz
popFile=/u/flashscratch/a/ab08028/captures/samples/samplesPop.Headers.forEasySFS.txt

easySFS=/u/home/a/ab08028/klohmueldata/annabel_data/bin/easySFS/easySFS.py

$easySFS -i $vcfdir/$vcf -p $popFile --preview


# test: needs to not be compressed?
gunzip snp_8_rmRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf.gz
