#! /bin/bash
#$ -cwd
#$ -l h_rt=3:00:00,h_data=16G,highp
#$ -N vcf1e_pullOutNeutral
#$ -o /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -e /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -m abe
#$ -M ab08028

# can also be run in the shell in about 1 hours + 30 mins
source /u/local/Modules/default/init/modules.sh

module load java
module load bedtools
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar

noCallFrac=1.0 # no filter used
maxHetFilter=0.75 # het filter used across all samples (per population het filter occurs during easy sfs)
#### file locations
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/vcf_filtering
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/sea_otter_genome/dedup_99_indexed_USETHIS/sea_otter_23May2016_bS9RH.deduped.99.fasta
REFSHORTCODE=SSO
#### parameters:
rundate=20200719_${REFSHORTCODE} # date genotypes were called and ref code 20200719_NSO 

vcfdir=$wd/${rundate}_filtered # date you called genotypes
outdir=$vcfdir/neutral_and_cds_VCFs/neutralVCFs
mkdir -p $outdir
 # note that these still have per-population all-het sites present ; but I modified the easy sfs script to remove them ()
# 20181018 bed file of neutral sites that have been called (min 10kb from genes, not in CpG island, doesn't blast to zebra fish); max no call frac is 0.9 (liberal)
neutralBed=/u/home/a/ab08028/klohmueldata/annabel_data/sea_otter_genome/sequenceCaptureTargetBedfiles/finalFiles-0based-Correction_20170906_USEFROMNOWON/finalNeutralSet.10K.1kb.regions.NoOverhang.NoFishMatch.NoNNNs.DupScaff99Removed.0-based.20170906.bed
# these are my neutral targets that I designed for SSO

# pull out neutral regions from the all_8 file (to use for getting monomorphic sites) and from the snp_8 file (for easy sfs projection)
#allVCF=all_8_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters_raw_variants.vcf.gz
allVCF=all_9_maxHetFilter_${maxHetFilter}_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters_raw_variants.vcf.gz

#snpVCF=snp_8b_forEasySFS_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf.gz
snpVCF=snp_9b_forEasySFS_maxHetFilter_${maxHetFilter}_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf.gz
# snp vcf first:  # for easy SFS, don't have snp sfs be gzipped 

java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/${snpVCF} \
-o $outdir/neutral.${snpVCF%.gz} \
-L $neutralBed

# then all-sfs  (DO output as .gz)
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/${allVCF} \
-o $outdir/neutral.${allVCF} \
-L $neutralBed

