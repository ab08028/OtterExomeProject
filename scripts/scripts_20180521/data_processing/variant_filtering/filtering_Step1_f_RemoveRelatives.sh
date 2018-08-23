#! /bin/bash
#$ -cwd
#$ -l h_rt=50:00:00,h_data=16G,highp
#$ -N vcf1f_removeRelatives
#$ -o /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -e /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -m abe
#$ -M ab08028
######## These indivduals were identified in PLINK as having a kinship coeff > 0.2 with another individual in the same population
# I am removing one from each pair 
# modules
source /u/local/Modules/default/init/modules.sh
module load java
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar

#### parameters:
rundate=20180806 # date genotypes were called (vcf_20180806 includes capture 02)
noCallFrac=0.2 # maximum fraction of genotypes that can be "no call" (./.) # note that genotypes can still "PASS" if 


#### file locations
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/vcf_filtering
mkdir -p $wd
vcfdir=$SCRATCH/captures/vcfs/vcf_${rundate}/
infile=raw_variants.vcf.gz ### make sure this doesn't have a path as part of its name! just infile names
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta

### Relatives to remove
ind1="RWAB003_19_ELUT_CA_352"
ind2="106_Elut_AL_AD_GE91109"
ind3="101_Elut_KUR_3"
ind4="102_Elut_KUR_4"
ind5="104_Elut_KUR_7"
ind6="77_Elut_KUR_14"

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
-xl_sn ${ind6}

# remove from snp file:
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/'snp_8_maxNoCallFrac_'${noCallFrac}'_passingBespoke_passingAllFilters_postMerge_'${infile} \
-o ${vcfdir}/'snp_8_rmRelativesAdmixed_maxNoCallFrac_'${noCallFrac}'_passingBespoke_passingAllFilters_postMerge_'${infile} \
-xl_sn ${ind1} \
-xl_sn ${ind2} \
-xl_sn ${ind3} \
-xl_sn ${ind4} \
-xl_sn ${ind5} \
-xl_sn ${ind6}

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
-xl_sn ${ind6}

############### remake population vcfs: ####################
# Bering/Medny --> Commanders (COM)
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${outdir}/'all_8_rmRelativesAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
-o ${outdir}/populationVCFs/COM_'all_8_rmRelativesAdmixed_passingAllFilters_allCalled.vcf.gz' \
-se '.+_Elut_BER_.+' \
-se '.+_Elut_MED_.+' \
--maxNOCALLfraction 0

# California --> CA (include the RWAB hiseq4000 samples)
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${outdir}/'all_8_rmRelativesAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
-o ${outdir}/populationVCFs/CA_'all_8_rmRelativesAdmixed_passingAllFilters_allCalled.vcf.gz' \
-se '.+_Elut_CA_.+' \
-se 'RWAB003_.+_ELUT_CA_.+' \
--maxNOCALLfraction 0

# Alaska --> AK
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${outdir}/'all_8_rmRelativesAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
-o ${outdir}/populationVCFs/AK_'all_8_rmRelativesAdmixed_passingAllFilters_allCalled.vcf.gz' \
-se '.+_Elut_AK_.+' \
--maxNOCALLfraction 0

# Aleutian --> AL
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${outdir}/'all_8_rmRelativesAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
-o ${outdir}/populationVCFs/AL_'all_8_rmRelativesAdmixed_passingAllFilters_allCalled.vcf.gz' \
-se '.+_Elut_AL_.+' \
--maxNOCALLfraction 0


# Kuril --> KUR (include the RWAB hiseq4000 samples)
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${outdir}/'all_8_rmRelativesAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
-o ${outdir}/populationVCFs/KUR_'all_8_rmRelativesAdmixed_passingAllFilters_allCalled.vcf.gz' \
-se '.+_Elut_KUR_.+' \
-se 'RWAB003_.+_ELUT_KUR_.+' \
--maxNOCALLfraction 0