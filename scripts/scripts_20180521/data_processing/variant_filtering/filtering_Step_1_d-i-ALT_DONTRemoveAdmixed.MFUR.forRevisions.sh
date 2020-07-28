#! /bin/bash
#$ -cwd
#$ -l h_rt=30:00:00,h_data=16G,highp
#$ -N vcf1d_removeRelatives
#$ -o /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -e /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -m abe
#$ -M ab08028

####### Need a version of filter 1di that doesn't remove admixed
source /u/local/Modules/default/init/modules.sh
module load java
module load python
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar
tabix=/u/home/a/ab08028/klohmueldata/annabel_data/bin/tabix-0.2.6/tabix
bgzip=/u/home/a/ab08028/klohmueldata/annabel_data/bin/tabix-0.2.6/bgzip

noCallFrac=0.2 # ### adding this in for per individual het ; this is new for the per ind het calculations. max missingness (./.) per genotype allowed (so 20% can be missing at a site but more than that means its filtered out)
#snpNoCallFrac=0.2 # max missingness allowed in snp file (stricter cutoff)
maxHetFilter=0.75 # (maximum heterozygous genotypes per site filter)

REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta
#### parameters:
rundate=20181119 # date genotypes were called (vcf_20180806 includes capture 02)

#### file locations
#SCRATCH=/u/flashscratch/a/ab08028
wd=/u/home/a/ab08028/klohmueldata/annabel_data/captures/vcf_filtering

mkdir -p $wd
infile=raw_variants.vcf.gz ### make sure this doesn't have a path as part of its name! just infile names

vcfdir=$wd/${rundate}_filtered # date you called genotypes

gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/
scriptdir=$gitdir/data_processing/variant_filtering
hetFilterScript=filtering_perPopulation.noCall.maxHetFilter.py
######### KEEPING THESE ALL IN in this alternate script for revisions #########
## Relatives to remove
ind1="RWAB003_19_ELUT_CA_352"
ind2="106_Elut_AL_AD_GE91109"
#ind3="101_Elut_KUR_3"
#ind4="102_Elut_KUR_4"
# switched ind3 and ind4 to remove other member of relative pair (removing RWAB) -- 20181204
#ind3="RWAB003_28_ELUT_KUR_12" # already removed RWAB samples
#ind4="RWAB003_26_ELUT_KUR_15"  # already removed RWAB samples
ind5="104_Elut_KUR_7"
ind6="77_Elut_KUR_14"

########## kEEP THESE IN:
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
ind19="141_Elut_CA_419" # doesn't appear admixed but is a PCA outlier 
###### ** UPDATED THIS LIST *** for 20181119 ## 

## 20180910 update: I am not going to remove admixed individuals from the all_8 vcf files
# Because they are valid sequences for PCA, etc. Just removing them from pop-specific files

############################# remove admixed and relatives and other outliers from the all file ###############
# update in 20181206: not making an all file that removes relatives but keeps admixed; that is excessive file size
# and you can get that from all_7 if you need it (just exclude relatives)
# java -jar $GATK \
# -R $REFERENCE \
# -T SelectVariants \
# --variant ${vcfdir}/'all_7_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
# -o ${vcfdir}/'all_8_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
# -xl_sn ${ind1} \
# -xl_sn ${ind2} \
# -xl_sn ${ind5} \
# -xl_sn ${ind6} \
# -xl_sn ${ind7} \
# -xl_sn ${ind8} \
# -xl_sn ${ind9} \
# -xl_sn ${ind10} \
# -xl_sn ${ind11} \
# -xl_sn ${ind12} \
# -xl_sn ${ind13} \
# -xl_sn ${ind14} \
# -xl_sn ${ind15} \
# -xl_sn ${ind16} \
# -xl_sn ${ind17} \
# -xl_sn ${ind18} \
# -xl_sn ${ind19} \
# -trimAlternates
# # adding trim alternates? 20181214


######## removing all outliers, relatives and admixed here ########

#######################################################################################
###################################### imposing excess heterozygosity filter ################################
#######################################################################################

##### going to filter out sites that exceed a number of het-calls across all samples (then I do more specific population filtering)
# of heterozygosity within easysfs. But the goal of this stage is to make sure that my files that go into PCA, etc. aren't too affected by these
# weird sites.
#######  this will be file for heterozygosity so want to add in a missingness filter 
# going straight from all _7 instead of all_8 bc i want the admixed in there
python $scriptdir/$hetFilterScript --vcf ${vcfdir}/'all_7_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters_'${infile} \
--outfile ${vcfdir}/'all_10_maxHetFilter_'${maxHetFilter}'_ForRevPerIndHet_keepRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile%.gz} \
--errorfile ${vcfdir}/'sitesFailingMaxHetFilter_'${maxHetFilter}'.txt' \
--maxNoCallFrac $noCallFrac \
--maxHetFilter $maxHetFilter

$bgzip -f ${vcfdir}/'all_10_maxHetFilter_'${maxHetFilter}'_ForRevPerIndHet_keepRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile%.gz}
$tabix -f -p vcf ${vcfdir}/'all_10_maxHetFilter_'${maxHetFilter}'_ForRevPerIndHet_keepRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile}


# then can pull neutral sites out of this. 


