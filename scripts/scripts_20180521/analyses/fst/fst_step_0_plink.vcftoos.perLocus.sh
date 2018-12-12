########### You know, I think it's best to do this in SNP Relate in R <-------- ##########

module load vcftools
module load plink
gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/
scriptdir=$gitdir/analyses/TREEMIX/
genotypeDate=20181119
vcfdir=/u/flashscratch/a/ab08028/captures/vcf_filtering/${genotypeDate}_filtered/

# eventually want to use the "clean" vcf file (but I think it doesn't matter a ton because I'm only selecting "clean" individuals out)
header=snp_7_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants
fstdir=/u/flashscratch/a/ab08028/captures/analyses/fst/$genotypeDate/


##### fst with plink ######
plink --bfile $plinkFileDir/$header \
--freq \
--missing \
--within $clusters \
--allow-extra-chr \
--out $fstdir/plink.$header \
--nonfounders \
--keep-allele-order \
--fst

##### fst with vcftools ######

# CA.txt is file with list of CA individuals in it.

vcftools --gzvcf $vcfdir/$header.vcf.gz --weir-fst-pop CA.txt --weir-fst-pop BAJ.txt --out $fstdir/vcftools_CA_vs_BAJ
gzip -f $fstdirvcftools_CA_vs_BAJ.txt
# CA vs AK
vcftools --gzvcf $vcfdir/$header.vcf.gz --weir-fst-pop CA.txt --weir-fst-pop AK.txt --out $fstdir/vcftools_CA_vs_AK
gzip -f $fstdirvcftools_CA_vs_AK.txt
# AK vs BAJ
vcftools --gzvcf $vcfdir/$header.vcf.gz --weir-fst-pop BAJ.txt --weir-fst-pop AK.txt --out $fstdir/vcftools_BAJ_vs_AK
gzip -f $fstdirvcftools_BAJ_vs_AK.txt

# AK vs. AL

# AK. vs. KUR

# AK vs COM

# AL vs. COM

# AL vs. AK

########### You know, I think it's best to do this in SNP Relate in R <-------- ##########


# from log files:
#  get average across loci (mean Fst)
# better practice to average HS across loci, then HT across loci, ( I believe this is the "weighted" Fst)
# vcftools gives this info in the .log file. 

