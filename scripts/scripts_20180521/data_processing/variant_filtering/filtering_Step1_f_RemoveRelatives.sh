#! /bin/bash
#$ -cwd
#$ -l h_rt=50:00:00,h_data=16G,highp
#$ -N vcf1f_removeRelatives
#$ -o /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -e /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -m abe
#$ -M ab08028
######## These indivduals were identified in PLINK as having a kinship coeff > 0.2 with another individual in the same population
### you need to manually enter which individuals to remove; this is not an automated script; so be careful!
# I am removing one from each pair 
# modules
source /u/local/Modules/default/init/modules.sh
module load java
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar

#### parameters:
rundate=20180806 # date genotypes were called (vcf_20180806 includes capture 02)
noCallFrac=0.9 # maximum fraction of genotypes that can be "no call" (./.) that was used in previous steps in early files
allSitesNoCallFrac=0.2 # max missingness allowed in final file for all pops together (20%)
perPopNoCallFrac=0 # max missingness allowed in final file for each pop for sfs

#### file locations
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/vcf_filtering
mkdir -p $wd
infile=raw_variants.vcf.gz ### make sure this doesn't have a path as part of its name! just infile names
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta
vcfdir=$wd/${rundate}_filtered # date you called genotypes
mkdir -p $wd/populationVCFs
mkdir -p $wd/populationVCFs/admixedVCFs

## Relatives to remove
ind1="RWAB003_19_ELUT_CA_352"
ind2="106_Elut_AL_AD_GE91109"
ind3="101_Elut_KUR_3"
ind4="102_Elut_KUR_4"
ind5="104_Elut_KUR_7"
ind6="77_Elut_KUR_14"

### Heavily admixed invididuals / PCA outliers to remove:
ind7="90_Elut_KUR_20"
ind8="80_Elut_KUR_19"
ind9="91_Elut_KUR_21"
ind10="92_Elut_KUR_22"
ind11="93_Elut_KUR_23"
ind12="68_Elut_AK_AF3370"
ind13="56_Elut_AK_AF24903"
ind14="57_Elut_AK_AF24906"
ind15="163_Elut_AK_AF24915"
## Heavily admixed invididuals / PCA outliers to remove:

# Samples to exclude:
#RWAB003_19_ELUT_CA_352 (was paired with RWAB003_20_ELUT_CA_355; kinship 0.216347723514731 )
#106_Elut_AL_AD_GE91109 ( was paired with 118_Elut_AL_AD_GE91101; kinship 0.4845482485185 )
#101_Elut_KUR_3 ( was paired with RWAB003_28_ELUT_KUR_12; kinship  0.479370763531055 )
#102_Elut_KUR_4 ( was paired with RWAB003_26_ELUT_KUR_15; kinship 0.468276048119186 )
#104_Elut_KUR_7 ( was paired with 79_Elut_KUR_17; kinship 0.48474476630868 )
#77_Elut_KUR_14 ( was paired with 81_Elut_KUR_2; kinship 0.483877991667127 )
# Removing them from several files:
# remove from all sites file:
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/'all_7_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
-o ${vcfdir}/'all_8_rmRelativesAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
-xl_sn ${ind1} \
-xl_sn ${ind2} \
-xl_sn ${ind3} \
-xl_sn ${ind4} \
-xl_sn ${ind5} \
-xl_sn ${ind6} \
-xl_sn ${ind7} \
-xl_sn ${ind8} \
-xl_sn ${ind9} \
-xl_sn ${ind10} \
-xl_sn ${ind11} \
-xl_sn ${ind12} \
-xl_sn ${ind13} \
-xl_sn ${ind14} \
-xl_sn ${ind15}

# use the "8" version of all-sites for making SFSs
# but also make  a "9" version where you have 20% missingness cut off?

# remove from snp file:
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/'snp_7_maxNoCallFrac_'${noCallFrac}'_passingBespoke_passingAllFilters_postMerge_'${infile} \
-o ${vcfdir}/'snp_8_rmRelativesAdmixed_maxNoCallFrac_'${noCallFrac}'_passingBespoke_passingAllFilters_postMerge_'${infile} \
-xl_sn ${ind1} \
-xl_sn ${ind2} \
-xl_sn ${ind3} \
-xl_sn ${ind4} \
-xl_sn ${ind5} \
-xl_sn ${ind6} \
-xl_sn ${ind7} \
-xl_sn ${ind8} \
-xl_sn ${ind9} \
-xl_sn ${ind10} \
-xl_sn ${ind11} \
-xl_sn ${ind12} \
-xl_sn ${ind13} \
-xl_sn ${ind14} \
-xl_sn ${ind15}

# remove from nv file:
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/'nv_7_maxNoCallFrac_'${noCallFrac}'_passingBespoke_passingAllFilters_postMerge_'${infile} \
-o ${vcfdir}/'snp_8_rmRelativesAdmixed_maxNoCallFrac_'${noCallFrac}'_passingBespoke_passingAllFilters_postMerge_'${infile} \
-xl_sn ${ind1} \
-xl_sn ${ind2} \
-xl_sn ${ind3} \
-xl_sn ${ind4} \
-xl_sn ${ind5} \
-xl_sn ${ind6} \
-xl_sn ${ind7} \
-xl_sn ${ind8} \
-xl_sn ${ind9} \
-xl_sn ${ind10} \
-xl_sn ${ind11} \
-xl_sn ${ind12} \
-xl_sn ${ind13} \
-xl_sn ${ind14} \
-xl_sn ${ind15}


############## also put admixed individuals (but not relatives) into their own VCFs for later ##############

## KURIL ADMIXED:

java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/'all_7_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
-o ${vcfdir}/populationVCFs/admixedVCFs/admixIndOnly_KUR_'all_8_rmRelatives_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_allCalled.vcf.gz' \
-sn ${ind7} \
-sn ${ind8} \
-sn ${ind9} \
-sn ${ind10} \
-sn ${ind11} \
--maxNOCALLfraction $perPopNoCallFrac

## Alaska ADMIXED:

java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/'all_7_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
-o ${vcfdir}/populationVCFs/admixedVCFs/admixIndOnly_AK_'all_8_rmRelatives_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_allCalled.vcf.gz' \
-sn ${ind12} \
-sn ${ind13} \
-sn ${ind14} \
-sn ${ind15} \
--maxNOCALLfraction $perPopNoCallFrac


############### remake population vcfs without relatives or admixed: ####################
# Bering/Medny --> Commanders (COM)
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/'all_8_rmRelativesAdmixed_passingBespoke_rmBadIndividuals_passingFilters_'${infile} \
-o ${vcfdir}/populationVCFs/COM_'all_8_rmRelativesAdmixed_passingAllFilters_allCalled.vcf.gz' \
-se '.+_Elut_BER_.+' \
-se '.+_Elut_MED_.+' \
--maxNOCALLfraction $perPopNoCallFrac

# California --> CA (include the RWAB hiseq4000 samples)
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/'all_8_rmRelativesAdmixed_passingBespoke_rmBadIndividuals_passingFilters_'${infile} \
-o ${vcfdir}/populationVCFs/CA_'all_8_rmRelativesAdmixed_passingAllFilters_allCalled.vcf.gz' \
-se '.+_Elut_CA_.+' \
-se 'RWAB003_.+_ELUT_CA_.+' \
--maxNOCALLfraction $perPopNoCallFrac

# Alaska --> AK
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/'all_8_rmRelativesAdmixed_passingBespoke_rmBadIndividuals_passingFilters_'${infile} \
-o ${vcfdir}/populationVCFs/AK_'all_8_rmRelativesAdmixed_passingAllFilters_allCalled.vcf.gz' \
-se '.+_Elut_AK_.+' \
--maxNOCALLfraction $perPopNoCallFrac

# Aleutian --> AL
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/'all_8_rmRelativesAdmixed_passingBespoke_rmBadIndividuals_passingFilters_'${infile} \
-o ${vcfdir}/populationVCFs/AL_'all_8_rmRelativesAdmixed_passingAllFilters_allCalled.vcf.gz' \
-se '.+_Elut_AL_.+' \
--maxNOCALLfraction $perPopNoCallFrac


# Kuril --> KUR (include the RWAB hiseq4000 samples)
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/'all_8_rmRelativesAdmixed_passingBespoke_rmBadIndividuals_passingFilters_'${infile} \
-o ${vcfdir}/populationVCFs/KUR_'all_8_rmRelativesAdmixed_passingAllFilters_allCalled.vcf.gz' \
-se '.+_Elut_KUR_.+' \
-se 'RWAB003_.+_ELUT_KUR_.+' \
--maxNOCALLfraction $perPopNoCallFrac

######################### after pop-vcfs made, make a final version of the all sites vcfs that uses a 20% cutoff -- use these for PCA, fastructure,etc. #########
# snps file: 
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/'snp_8_rmRelativesAdmixed_maxNoCallFrac_'${noCallFrac}'_passingBespoke_passingAllFilters_postMerge_'${infile}
-o ${vcfdir}/'snp_9_rmRelativesAdmixed_maxNoCallFrac_'${allSitesNoCallFrac}'_passingBespoke_passingAllFilters_postMerge_'${infile} \
--maxNOCALLfraction ${allSitesNoCallFrac}

# nv file: 
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/'nv_8_rmRelativesAdmixed_maxNoCallFrac_'${noCallFrac}'_passingBespoke_passingAllFilters_postMerge_'${infile}
-o ${vcfdir}/'nv_9_rmRelativesAdmixed_maxNoCallFrac_'${allSitesNoCallFrac}'_passingBespoke_passingAllFilters_postMerge_'${infile} \
--maxNOCALLfraction ${allSitesNoCallFrac}


# all sites file: 
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/'all_8_rmRelativesAdmixed_maxNoCallFrac_'${noCallFrac}'_passingBespoke_passingAllFilters_postMerge_'${infile}
-o ${vcfdir}/'all_9_rmRelativesAdmixed_maxNoCallFrac_'${allSitesNoCallFrac}'_passingBespoke_passingAllFilters_postMerge_'${infile} \
--maxNOCALLfraction ${allSitesNoCallFrac}

