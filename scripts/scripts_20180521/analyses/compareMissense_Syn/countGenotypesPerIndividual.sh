source /u/local/Modules/default/init/modules.sh
module load python/2.7
# Count up types of sites for missense and synonymous and neutral

genotypeDate=20181119
wd=$SCRATCH/captures/analyses/compareMissense_Syn
indir=/u/flashscratch/a/ab08028/captures/vcf_filtering/${genotypeDate}_filtered/neutral_and_cds_VCFs/
gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject
scriptdir=$gitdir/scripts/scripts_20180521/analyses/compareMissense_Syn
script=countHomRefHomAltHetNoCall.perIndividual.py
todaysdate=`date +%Y%m%d`
############################ cds -- annotated SNPS only! ####################
# vep was carried out during filtering steps (1e)
# so now have missense and synonymous separate vcfs 

outdir=$wd/countsOfGenotypesPerIndividual
mkdir -p $outdir
# can count up the number of hets/homs per individual for each one
vcfIdentifier=snp_9b_forEasySFS_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf
echo "These are the vcfs being used to get the counts:" > $outdir/countsLog.vcfsUsed.${todaysdate}.txt

# note full vcf name is missense_vep_cds_snp_9b_forEasySFS_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf
for category in missense syn
do
echo "starting $category"
vcf=${category}_vep_cds_${vcfIdentifier}
echo "$category: $vcf" >> $outdir/countsLog.vcfsUsed.${todaysdate}.txt
python $scriptdir/$script --vcf $indir/cdsVCFs/$vcf --outdir $outdir --outPREFIX ${category}.countsPerIndividual
done


################# neutral is different -- all sites not just SNPs (to get heterozygosity) #######################
# for getting heterozygosity you need to use the full vcf file 

category=neutral
vcfIdentifier=all_9_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters_raw_variants.vcf.gz
vcf=${category}.${vcfIdentifier}
echo "starting neutral"
echo "$category: $vcf" >> $outdir/countsLog.vcfsUsed.${todaysdate}.txt
# is in a different directory, but output goes to same directory
python $scriptdir/$script --vcf $indir/neutralVCFs/$vcf --outdir $outdir --outPREFIX ${category}.countsPerIndividual
