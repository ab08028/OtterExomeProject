#! /bin/bash
#$ -cwd
#$ -l h_rt=50:00:00,h_data=8G,highp
#$ -N countSitesNEW
#$ -o /u/flashscratch/a/ab08028/captures/reports/cdsSites
#$ -e /u/flashscratch/a/ab08028/captures/reports/cdsSites
#$ -m abe
#$ -M ab08028
##### testing: see how missing data factors in
# keep only sites where all are called (and then also maybe try 80% called)

source /u/local/Modules/default/init/modules.sh
module load python/2.7
module load vcftools
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


# can count up the number of hets/homs per individual for each one
vcfIdentifier=snp_9b_forEasySFS_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_passingBespoke_passingAllFilters_postMerge_raw_variants
echo "These are the vcfs being used to get the counts:" > $outdir/countsLog.vcfsUsed.${todaysdate}.txt

# note full vcf name is missense_vep_cds_snp_9b_forEasySFS_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf
for minCallRate in 1 0.8 # 1 = no missing allowed; 0.8 = 20% missing allowed 
do
outdir=$wd/countsOfGenotypesPerIndividual/minCallRate_$minCallRate
mkdir -p $outdir
echo "starting min call rate $minCallRate"
for category in missense syn
do
echo "starting $category"
vcf=${category}_vep_cds_${vcfIdentifier}.vcf.gz

# remove all sites with any missing data (should drop sites down a lot)
# note weird naming of max-missing in vcftools; it actually refers to the min call rate allowed so 0 is all sites allowed
# 1 is no missing data allowed; 0.8 says that at least 0.8 must be called 
vcftools --gzvcf $indir/cdsVCFs/$vcf \
	--max-missing $minCallRate \
	--out $indir/cdsVCFs/${category}_vep_cds_${vcfIdentifier}.minCall.$minCallRate \
	--recode
noMissingVCF=${category}_vep_cds_${vcfIdentifier}.minCall.$minCallRate.recode.vcf
echo "$category: $noMissingVCF" >> $outdir/countsLog.vcfsUsed.${todaysdate}.txt
python $scriptdir/$script --vcf $indir/cdsVCFs/$noMissingVCF --outdir $outdir --outPREFIX ${category}.countsPerIndividual.minCall.$minCallRate
done
done
# then want to do the counts:

# maybe want to remove totally monomorphic sites too? try that next perhaps

############### count total cds sites: #########
for minCallRate in 1 0.8
do
echo "starting rate $minCallRate"
vcf=cds_all_9_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters_raw_variants.vcf.gz
#vcftools --gzvcf $indir/cdsVCFs/$vcf \
#--max-missing $minCallRate 	\
#--out $indir/cdsVCFs/${vcf%.vcf.gz}.minCallRate.${minCallRate} \
#--recode

gzip -f $indir/cdsVCFs/${vcf%.vcf.gz}.minCallRate.${minCallRate}.recode.vcf
# count callable sites per ind:
script=countCallableSitesPerIndividual.py
outdir=/u/flashscratch/a/ab08028/captures/analyses/compareMissense_Syn/callableCDSPerIndividual/minCallRate_$minCallRate
mkdir -p $outdir
### run script: 
python $scriptdir/$script --vcf $indir/cdsVCFs/${vcf%.vcf.gz}.minCallRate.${minCallRate}.recode.vcf.gz --outfile $outdir/callableCDSSitesPerIndividual.${todaysdate}.txt

done
# this may need to be a job (make take a while)