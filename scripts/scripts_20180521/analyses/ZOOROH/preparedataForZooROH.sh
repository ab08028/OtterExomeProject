############## Want to try ZooROH #############
####
# need to convert to "Oxford Gen" format (pg 10 of manual)
module load bcftools
vcf=snp_9b_forEasySFS_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf.gz

bcftools convert -t ^chrX,chrY,chrM -g outfile --chrom --tag GT myfile.vcf

# convert missing GTs into 0
sed -e 's/-nan/0/g' file.gen > newfile.gen
