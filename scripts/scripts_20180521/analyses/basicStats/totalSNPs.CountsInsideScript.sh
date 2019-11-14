
#### count total snps, excluding monomorphic

source /u/local/Modules/default/init/modules.sh
module load python/2.7
module load vcftools
genotypeDate=20181119
wd=$SCRATCH/captures/analyses/compareMissense_Syn
indir=/u/flashscratch/a/ab08028/captures/vcf_filtering/${genotypeDate}_filtered//
gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject
scriptdir=$gitdir/scripts/scripts_20180521/analyses/compareMissense_Syn
script=countHomRefHomAltHetNoCall.perIndividual.py
todaysdate=`date +%Y%m%d`

# use all_7 because want to get GT counts for admixed individuals too
# want to count snps with and without 1/1 sites
# wait I don't need all then do I.
# just need to count snp files 
# only need all for total called genotypes per individual 

allVCF=all_8_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters_raw_variants.vcf.gz
# use this for total GT counts per individual


# but for counting snps use snp_9a 
# this is the number of snps with a 0.2 missingness filter, excluding sites that are monomorphic across all inds
# 89 individuals

vcftools \
	--gzvcf $indir/snp_9a_forPCAetc_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf.gz \
	--mac 1 \
	--out /u/flashscratch/a/ab08028/captures/basicSNPstats/countingSNPs.fromSNP9a.${todaysdate}.vcftools.out

vcftools \
	--gzvcf $indir/snp_7_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf.gz \
	--mac 1 \
	--out /u/flashscratch/a/ab08028/captures/basicSNPstats/countingSNPs.fromSNP7.${todaysdate}.vcftools.out
# need to do for snp_7 as well to include admixed 
# have the log go there ^
# result (2019114)

#After filtering, kept 89 out of 89 Individuals
#After filtering, kept 52628 out of a possible 1571275 Sites
#Run Time = 88.00 seconds
# so it's 1571275 sites when including 1/1 sites
# and it's only 52628 sites that are polymorphic among otters (all pops)

# checked it using grep:
#zcat snp_9a_forPCAetc_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf.gz | grep -v "#" | grep "1/1" | grep -v "0/0" | grep -v "0/1" -c
#yields 1518647 sites that are entirely 1/1 not with any 0/0 or 0/1s
#this is eactly 1571275 (total) - 52628, so we are good to go on the snp count
# want to get counts for admixed (all_7) but also with all the final filters (all_9)
# shouldn't be that different, could probably plot either
# do I also want it for all_5 which shows why I removed bad individuals I guess
# sort of tied to the missingness plot I guess... 
allVCFs="all_7_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters_raw_variants.vcf.gz all_9_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters_raw_variants.vcf.gz"

# just want to report the number of called genotypes per individual after filtering (but do want admixed to be included in that?)
outdir=/u/flashscratch/a/ab08028/captures/basicSNPstats
for allVCF in $allVCFs
do
# and then also want to count total called genotypes per individual :
python $scriptdir/$script --vcf $indir/$allVCF --outdir $outdir --outPREFIX ${allVCF%.vcf.gz}.totalCalledGenotypesPerIndividual.${todaysdate}.txt
done
