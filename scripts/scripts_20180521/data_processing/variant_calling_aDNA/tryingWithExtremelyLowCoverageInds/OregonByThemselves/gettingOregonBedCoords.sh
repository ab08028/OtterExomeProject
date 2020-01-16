######## working with mfur mapped only ######
module load vcftools
module load bedtools
############## get snp9b variable sites: ############

input=/u/flashscratch/a/ab08028/captures/vcf_filtering/20181119_filtered/snp_9b_forEasySFS_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf.gz
outdir=/u/flashscratch/a/ab08028/captures/vcf_filtering/20181119_filtered/variableSitesBedCoords
intersectdir=/u/flashscratch/a/ab08028/captures/aDNA-ModernComparison/angsd-pseudoHaps/20191216-ORs-pseudohaps/intersectWithModernData
alleleFreqdir=/u/flashscratch/a/ab08028/captures/analyses/alleleFreqsInCA_AK
mkdir -p $outdir
bedhead="#chrom\tstart0based\tend\tmarkerID\tempty5\tempty6\tempty7\tempty8\tempty9\tempty10\tempty11\tempty12"
haphead="chr\tpos\tmajor\tind0\tind1"
vcfhead=`zcat $input | grep -v "##"| grep -m1 "#"`
comboheaderHAP=`echo -e "$bedhead\t$haphead"`
comboheaderVCF=`echo -e "$bedhead\t$vcfhead"`
# get mac >=1 sites
vcftools --gzvcf $input --mac 1 --recode --out $outdir/snp_9b_forEasySFS_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_passingBespoke_passingAllFilters_postMerge_raw_variants.NOFIXEDALT # sites that have at least one minor allele

gzip $outdir/snp_9b_forEasySFS_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_passingBespoke_passingAllFilters_postMerge_raw_variants.NOFIXEDALT.recode.vcf

# get bed coords from that (with everything else as info )
echo -e "$comboheaderVCF" >  $outdir/snp_9b_forEasySFS_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_passingBespoke_passingAllFilters_postMerge_raw_variants.NOFIXEDALT.superfile.bed
zcat $outdir/snp_9b_forEasySFS_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_passingBespoke_passingAllFilters_postMerge_raw_variants.NOFIXEDALT.recode.vcf.gz | grep -v "#" | awk '{OFS="\t";print $1,$2-1,$2,$1"_"$2,".",".",".",".",".",".",".",".",$0}' | sed 's/\t\t/\t/g' | sed 's/\t$//g' >> $outdir/snp_9b_forEasySFS_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_passingBespoke_passingAllFilters_postMerge_raw_variants.NOFIXEDALT.superfile.bed

gzip $outdir/snp_9b_forEasySFS_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_passingBespoke_passingAllFilters_postMerge_raw_variants.NOFIXEDALT.superfile.bed

############ then want to convert oregon haplo file to superfile: ##################
hapdir=/u/flashscratch/a/ab08028/captures/aDNA-ModernComparison/angsd-pseudoHaps/20191216-ORs-pseudohaps
echo -e "$comboheaderHAP" >  $hapdir/oregon.coords.superfile.bed 
zcat $hapdir/angsdOut.mappedTomfur.haplo.gz | grep -v "chr" | awk '{OFS="\t";print $1,$2-1,$2,$1"_"$2,".",".",".",".",".",".",".",".",$0}' | sed 's/\t\t/\t/g' | sed 's/\t$//g'>> $hapdir/oregon.coords.superfile.bed 

################ intersect beds ##############
bedtools intersect -header -a $hapdir/oregon.coords.superfile.bed -b $outdir/snp_9b_forEasySFS_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_passingBespoke_passingAllFilters_postMerge_raw_variants.NOFIXEDALT.superfile.bed.gz -wa > $intersectdir/Oregon.Modern.mfur.intersection.superfile.bed


############# then want to get frequency of those alleles in CA and AK from modern data ##########
# CA:
vcftools --gzvcf $outdir/snp_9b_forEasySFS_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_passingBespoke_passingAllFilters_postMerge_raw_variants.NOFIXEDALT.recode.vcf.gz \
--keep /u/flashscratch/a/ab08028/captures/samples/keep.CA.20181119.txt --freq --out $alleleFreqdir/CA.AlleleFreqs --bed $intersectdir/Oregon.Modern.mfur.intersection.superfile.bed

# AK:
vcftools --gzvcf $outdir/snp_9b_forEasySFS_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_passingBespoke_passingAllFilters_postMerge_raw_variants.NOFIXEDALT.recode.vcf.gz \
--keep /u/flashscratch/a/ab08028/captures/samples/keep.AK.20181119.txt --freq --out $alleleFreqdir/AK.AlleleFreqs --bed $intersectdir/Oregon.Modern.mfur.intersection.superfile.bed
