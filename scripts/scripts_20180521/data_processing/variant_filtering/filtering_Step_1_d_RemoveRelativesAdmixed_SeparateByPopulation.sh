#! /bin/bash
#$ -cwd
#$ -l h_rt=50:00:00,h_data=16G,highp
#$ -N vcf1d_removeRelatives
#$ -o /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -e /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -m abe
#$ -M ab08028
# Realized that maxnocallfrac doesn't act on the individuals you've selected, but rather the whole original VCF!
# so was being far too stringent! 
# So what I do now is I make a population vcf without a nocallfraction set - it will just have the baseline maxnocallfrac 0.9
# that was in the original file it is being drawn from. This might be useful for projecting down in dadi.
# THEN I am going to filter the allCalled
######## These indivduals were identified in PLINK as having a kinship coeff > 0.2 with another individual in the same population
### you need to manually enter which individuals to remove; this is not an automated script; so be careful!
# I am removing one from each pair 
# modules
source /u/local/Modules/default/init/modules.sh
module load java
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar

#### parameters:
rundate=20180806 # date genotypes were called (vcf_20180806 includes capture 02)
noCallFrac=0.9 # maximum fraction of genotypes that can be "no call" (./.) that was used in previous steps in previous "all sites" files (lenient cutoff)
snpNoCallFrac=0.2 # max missingness allowed in snp file (stricter cutoff)
perPopNoCallFrac=0 # max missingness allowed in final file for each pop for sfs (super strict cutoff)

#### file locations
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/vcf_filtering
mkdir -p $wd
infile=raw_variants.vcf.gz ### make sure this doesn't have a path as part of its name! just infile names
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta
vcfdir=$wd/${rundate}_filtered # date you called genotypes
mkdir -p $vcfdir/populationVCFs
mkdir -p $vcfdir/populationVCFs/admixedVCFs

## Relatives to remove
ind1="RWAB003_19_ELUT_CA_352"
ind2="106_Elut_AL_AD_GE91109"
ind3="101_Elut_KUR_3"
ind4="102_Elut_KUR_4"
ind5="104_Elut_KUR_7"
ind6="77_Elut_KUR_14"

### Heavily admixed invididuals / PCA outliers to keep out of population vcfs (can stay in all inds vcf):
ind7="90_Elut_KUR_20"
ind8="80_Elut_KUR_19"
ind9="91_Elut_KUR_21"
ind10="92_Elut_KUR_22"
ind11="93_Elut_KUR_23"
ind12="68_Elut_AK_AF3370"
ind13="56_Elut_AK_AF24903"
ind14="57_Elut_AK_AF24906"
ind15="163_Elut_AK_AF24915"
## mildly admixed invididuals / PCA outliers to remove:
ind16="125_Elut_AK_AF3369" # (adding)
ind17="65_Elut_AK_GE91060" # (adding)
ind18="165_Elut_AK_AL4661" # (doesn't appear admixed but is PCA outlier; so am keeping out of pop file and out of admixed file)


## 20180910 update: I am not going to remove admixed individuals from the all_8 vcf files
# Because they are valid sequences for PCA, etc. Just removing them from pop-specific files


java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/'all_7_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
-o ${vcfdir}/'all_8_rmRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
-xl_sn ${ind1} \
-xl_sn ${ind2} \
-xl_sn ${ind3} \
-xl_sn ${ind4} \
-xl_sn ${ind5} \
-xl_sn ${ind6} 
# not removing admixed at this stage 
#######################################################################################
############### make population vcfs without relatives or admixed: ####################
#######################################################################################
# Bering/Medny --> Commanders (COM)
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/'all_8_rmRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
-o ${vcfdir}/populationVCFs/COM_'all_8_rmRelativesAdmixed_passingAllFilters_maxNoCallFrac_'${noCallFrac}'.vcf.gz' \
-se '.+_Elut_BER_.+' \
-se '.+_Elut_MED_.+'
# note: no more maxnocall frac at this stage; file will be big -- but in all_9 my script filters it anyway

# note I'm keeping the old final name so that it flows into other scripts easily.

# California --> CA (include the RWAB hiseq4000 samples)
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/'all_8_rmRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
-o ${vcfdir}/populationVCFs/CA_'all_8_rmRelativesAdmixed_passingAllFilters_maxNoCallFrac_'${noCallFrac}'.vcf.gz' \
-se '.+_Elut_CA_.+' \
-se 'RWAB003_.+_ELUT_CA_.+' 
# relative was removed in main file above, so don't need to remove it here. 
# Alaska --> AK : remove admixed and make pop specific vcf
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/'all_8_rmRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
-o ${vcfdir}/populationVCFs/AK_'all_8_rmRelativesAdmixed_passingAllFilters_maxNoCallFrac_'${noCallFrac}'.vcf.gz' \
-se '.+_Elut_AK_.+' \
-xl_sn ${ind12} \
-xl_sn ${ind13} \
-xl_sn ${ind14} \
-xl_sn ${ind15} \
-xl_sn ${ind16} \
-xl_sn ${ind17} \
-xl_sn ${ind18}

# Aleutian --> AL
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/'all_8_rmRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
-o ${vcfdir}/populationVCFs/AL_'all_8_rmRelativesAdmixed_passingAllFilters_maxNoCallFrac_'${noCallFrac}'.vcf.gz' \
-se '.+_Elut_AL_.+'



# Kuril --> KUR (include the RWAB hiseq4000 samples)
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/'all_8_rmRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
-o ${vcfdir}/populationVCFs/KUR_'all_8_rmRelativesAdmixed_passingAllFilters_maxNoCallFrac_'${noCallFrac}'.vcf.gz' \
-se '.+_Elut_KUR_.+' \
-se 'RWAB003_.+_ELUT_KUR_.+' \
-xl_sn ${ind7} \
-xl_sn ${ind8} \
-xl_sn ${ind9} \
-xl_sn ${ind10} \
-xl_sn ${ind11} 

#######################################################################################
############## also put admixed individuals (but not relatives) into their own VCFs for later ##############
#######################################################################################
## KURIL ADMIXED (you take out of the all_7 file that still has admixed inds. in there and you positively select them with sn):

java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/'all_8_rmRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
-o ${vcfdir}/populationVCFs/admixedVCFs/admixIndOnly_KUR_'all_8_passingAllFilters_maxNoCallFrac_'${noCallFrac}'.vcf.gz' \
-sn ${ind7} \
-sn ${ind8} \
-sn ${ind9} \
-sn ${ind10} \
-sn ${ind11} 
## Alaska ADMIXED:

java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/'all_8_rmRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
-o ${vcfdir}/populationVCFs/admixedVCFs/admixIndOnly_AK_'all_8_passingAllFilters_maxNoCallFrac_'${noCallFrac}'.vcf.gz' \
-sn ${ind12} \
-sn ${ind13} \
-sn ${ind14} \
-sn ${ind15} \
-sn ${ind16} \
-sn ${ind17} 
# skipping ind 18 because it odesn't appear admixed in FASTRUCTURE; just appears like a PCA outlier

#######################################################################################
######################### after pop-vcfs made, make a final version  #########################
######################### of the all sites vcfs that uses a 20% cutoff #######################
#######################################################################################
# snps file: 
# Some sites that may have started as variant may have become INVARIANT by the end.
# Want these to end up in the INVARIANT category. 
# select the biallelic snps: 
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-V ${vcfdir}/'all_8_rmRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
--restrictAllelesTo BIALLELIC \
--selectTypeToInclude SNP \
-o ${vcfdir}/'snp_8a_rmRelatives_keepAdmixedOutliers_passingBespoke_maxNoCallFrac_'${snpNoCallFrac}'_passingBespoke_passingAllFilters_postMerge_'${infile} \
--maxNOCALLfraction ${snpNoCallFrac} 
# this call to maxnocallfrac is okay because it's for the snp file with a 20% cutoff for use in pca, etc.
########### this also keeps PCA outliers.

#### new: 20181130 Also want a version that excludes admixed and outliers and doesn't impose 0.2 filter:  #######
# use for easy SFS, etc.
######## all sites #########
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-V ${vcfdir}/'all_8_rmRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
-o ${vcfdir}/'all_8a_rmRelatives_rmAdmixedOutliers_passingBespoke_maxNoCallFrac_'${noCallFrac}'_passingBespoke_passingAllFilters_postMerge_'${infile} \
--maxNOCALLfraction ${noCallFrac} \
-xl_sn ${ind7} \
-xl_sn ${ind8} \
-xl_sn ${ind9} \
-xl_sn ${ind10} \
-xl_sn ${ind11} \
-xl_sn ${ind12} \
-xl_sn ${ind13} \
-xl_sn ${ind14} \
-xl_sn ${ind15} \
-xl_sn ${ind16} \
-xl_sn ${ind17} \
-xl_sn ${ind18}


#### only snps #############
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-V ${vcfdir}/'all_8_rmRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
--restrictAllelesTo BIALLELIC \
--selectTypeToInclude SNP \
-o ${vcfdir}/'snp_8a_rmRelatives_rmAdmixedOutliers_passingBespoke_maxNoCallFrac_'${noCallFrac}'_passingBespoke_passingAllFilters_postMerge_'${infile} \
--maxNOCALLfraction ${noCallFrac} \
-xl_sn ${ind7} \
-xl_sn ${ind8} \
-xl_sn ${ind9} \
-xl_sn ${ind10} \
-xl_sn ${ind11} \
-xl_sn ${ind12} \
-xl_sn ${ind13} \
-xl_sn ${ind14} \
-xl_sn ${ind15} \
-xl_sn ${ind16} \
-xl_sn ${ind17} \
-xl_sn ${ind18}

