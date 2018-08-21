#! /bin/bash
#$ -cwd
#$ -l h_rt=50:00:00,h_data=16G,highp,arch=intel*
#$ -N vcf1d_filtering
#$ -o /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -e /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -m abe
#$ -M ab08028


# modules
source /u/local/Modules/default/init/modules.sh
module load java
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar

rundate=20180806
noCallFrac=0.2
#### file locations
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/vcf_filtering
infile=raw_variants.vcf.gz ### make sure this doesn't have a path as part of its name! just infile names
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta

outdir=$wd/${rundate}_filtered # date you called genotypes
#################################################################################
############################ Separate by population #############################
#################################################################################
# want to separate vcfs by population
mkdir -p $outdir/populationVCFs
populations="CA AL KUR AK"

# treat commanders specially

# need to combine Bering and Medny separately (and label as COM (commanders))
########### bering and medny: ###########
# note: I'm making it so that every individual in the population has to be called at the site,
# but across all the populations only 80% have to be called at the site. I think this makes sense

# Bering/Medny --> Commanders (COM)
#java -jar $GATK \
#-R $REFERENCE \
#-T SelectVariants \
#--variant ${outdir}/'all_7_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
#-o ${outdir}/populationVCFs/COM_'all_7_passingAllFilters_allCalled'${infile} \
#-se '.+_Elut_BER_.+' \
#-se '.+_Elut_MED_.+' \
#--maxNOCALLfraction 0

# California --> CA (include the RWAB hiseq4000 samples)
#java -jar $GATK \
#-R $REFERENCE \
#-T SelectVariants \
#--variant ${outdir}/'all_7_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
#-o ${outdir}/populationVCFs/CA_'all_7_passingAllFilters_allCalled'${infile} \
#-se '.+_Elut_CA_.+' \
#-se 'RWAB003_.+_ELUT_CA_.+' \
#--maxNOCALLfraction 0

# Alaska --> AK
#java -jar $GATK \
#-R $REFERENCE \
#-T SelectVariants \
#--variant ${outdir}/'all_7_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
#-o ${outdir}/populationVCFs/AK_'all_7_passingAllFilters_allCalled'${infile} \
#-se '.+_Elut_AK_.+' \
#--maxNOCALLfraction 0

# Aleutian --> AL
#java -jar $GATK \
#-R $REFERENCE \
#-T SelectVariants \
#--variant ${outdir}/'all_7_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
#-o ${outdir}/populationVCFs/AL_'all_7_passingAllFilters_allCalled'${infile} \
#-se '.+_Elut_AL_.+' \
#--maxNOCALLfraction 0


# Kuril --> KUR (include the RWAB hiseq4000 samples)
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${outdir}/'all_7_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
-o ${outdir}/populationVCFs/KUR_'all_7_passingAllFilters_allCalled'${infile} \
-se '.+_Elut_KUR_.+' \
-se 'RWAB003_.+_ELUT_KUR_.+' \
--maxNOCALLfraction 0

# using $pop inside regex doesn't work; going pop by pop 
