#! /bin/bash
#$ -cwd
#$ -l h_rt=50:00:00,h_data=16G,highp
#$ -N makeSFS_eachStage
#$ -o /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -e /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -m abe
#$ -M ab08028
############### --> ! need to figure out the all/no call at each stage bc maxnocallfrac doesn't work the way you thought!

# sep by population and only select neutral sites
source /u/local/Modules/default/init/modules.sh
module load java
module load python/2.7

GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar
#### parameters:
rundate=20180806 # date genotypes were called (vcf_20180806 includes capture 02)
noCallFrac=0.9 # maximum fraction of genotypes that can be "no call" (./.) that was used in previous steps in previous "all sites" files (lenient cutoff)
snpNoCallFrac=0.2 # max missingness allowed in snp file (stricter cutoff)
perPopNoCallFrac=0 # max missingness allowed in final file for each pop for sfs (super strict cutoff)
#### file locations
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/vcf_filtering
#wd=/u/home/a/ab08028/nobackup-kirk/annabel/captures/vcf_filtering/
mkdir -p $wd
infile=raw_variants.vcf.gz ### make sure this doesn't have a path as part of its name! just infile names
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta
vcfdir=$wd/${rundate}_filtered # date you called genotypes
neutralBed=${vcfdir}/bedCoords/all_8_rmRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_0.9_rmBadIndividuals_passingFilters.min10kb.fromExon.noCpGIsland.noRepeat.noFish.0based.sorted.merged.useThis.bed
######## sfs specific ########

gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptdir=$gitdir/scripts/scripts_20180521/analyses/generate_sfs/
script=generate.1D.SFS.FilterNoCallSimultaneously.py # this script will exclude any line that has no-call genotypes
# is different from regular script that will error out if there is a nocall genotype present because it indicates
# a problem with filtering


## Relatives to remove

### Heavily admixed invididuals / PCA outliers to keep out of population vcfs (can stay in all inds vcf):

##############################################
# don't use maxnocallfrac! Instead use sfs filtering script

mkdir -p ${vcfdir}/step-by-step-vcf-sfs
infile=raw_variants.vcf.gz
# eventually
for vcfFile in 'all_1_TrimAlt_'${infile} 'all_5_passingFilters_'${infile} 'all_6_rmBadIndividuals_passingFilters_'${infile} 'all_7_passingBespoke_maxNoCallFrac_0.9_rmBadIndividuals_passingFilters_'${infile} 'all_8_rmRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_0.9_rmBadIndividuals_passingFilters_'${infile} 'snp_2_Filter_TrimAlt_raw_variants.vcf.gz snp_3_Flagged_GQ_DP_GaTKHF_cluster_'${infile} 'snp_4_Filtered_GQ_DP_GaTKHF_cluster_'${infile} 'snp_8_rmRelativesAdmixed_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge'${infile}


do
# get stats?
echo "starting on " $vcfFile

# Make Commanders SFS
echo "COM"
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant $vcfdir/$vcfFile \
-o ${vcfdir}/step-by-step-vcf-sfs/COM.allCalled.neutralOnly.$vcfFile \
-se '.+_Elut_BER_.+' \
-se '.+_Elut_MED_.+' \
-L $neutralBed \
-xl_sn 161_Elut_MED_19 \
-xl_sn 76_Elut_BER_31 \
-xl_sn 86_Elut_MED_25 \
-xl_sn 131_Elut_BER_34


python $scriptdir/$script --vcf ${vcfdir}/step-by-step-vcf-sfs/COM.allCalled.neutralOnly.$vcfFile --pop COM --outdir ${vcfdir}/step-by-step-vcf-sfs/ --outPREFIX ${vcfFile%_${infile}}

# California --> CA (include the RWAB hiseq4000 samples)
echo "CA"
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant $vcfdir/$vcfFile \
-o ${vcfdir}/step-by-step-vcf-sfs/CA.allCalled.neutralOnly.$vcfFile \
-se '.+_Elut_CA_.+' \
-se 'RWAB003_.+_ELUT_CA_.+' \
-L $neutralBed \
-xl_sn RWAB003_15_ELUT_CA_159 \
-xl_sn RWAB003_22_ELUT_CA_410 \
-xl_sn RWAB003_19_ELUT_CA_352


python $scriptdir/$script --vcf ${vcfdir}/step-by-step-vcf-sfs/CA.allCalled.neutralOnly.$vcfFile --pop CA --outdir ${vcfdir}/step-by-step-vcf-sfs/ --outPREFIX ${vcfFile%_${infile}}


# Alaska --> AK : remove admixed and make pop specific vcf
echo "AK"
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant $vcfdir/$vcfFile \
-o ${vcfdir}/step-by-step-vcf-sfs/AK.allCalled.neutralOnly.$vcfFile \
-se '.+_Elut_AK_.+' \
-L $neutralBed \
-xl_sn 128_Elut_AK_AF3630 \
-xl_sn 150_Elut_AK_AF3631 \
-xl_sn 151_Elut_AK_AF3714 \
-xl_sn 152_Elut_AK_AF3720 \
-xl_sn 153_Elut_AK_AF3181 \
-xl_sn 68_Elut_AK_AF3370 \
-xl_sn 56_Elut_AK_AF24903 \
-xl_sn 57_Elut_AK_AF24906 \
-xl_sn 163_Elut_AK_AF24915 \
-xl_sn 125_Elut_AK_AF3369 \
-xl_sn 65_Elut_AK_GE91060 \
-xl_sn 165_Elut_AK_AL4661


python $scriptdir/$script --vcf ${vcfdir}/step-by-step-vcf-sfs/AK.allCalled.neutralOnly.$vcfFile --pop AK --outdir ${vcfdir}/step-by-step-vcf-sfs/ --outPREFIX ${vcfFile%_${infile}}

# Aleutian --> AL
echo "AL"
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant $vcfdir/$vcfFile \
-o ${vcfdir}/step-by-step-vcf-sfs/AL.allCalled.neutralOnly.$vcfFile \
-se '.+_Elut_AL_.+' \
-L $neutralBed \
-xl_sn 112_Elut_AL_AM_AM92022 \
-xl_sn 106_Elut_AL_AD_GE91109

python $scriptdir/$script --vcf ${vcfdir}/step-by-step-vcf-sfs/AL.allCalled.neutralOnly.$vcfFile --pop AL --outdir ${vcfdir}/step-by-step-vcf-sfs/ --outPREFIX ${vcfFile%_${infile}}

# Kuril --> KUR (include the RWAB hiseq4000 samples)
echo "KUR"
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant $vcfdir/$vcfFile \
-o ${vcfdir}/step-by-step-vcf-sfs/KUR.allCalled.neutralOnly.$vcfFile \
-se '.+_Elut_KUR_.+' \
-se 'RWAB003_.+_ELUT_KUR_.+' \
-L $neutralBed \
-xl_sn 101_Elut_KUR_3 \
-xl_sn 102_Elut_KUR_4 \
-xl_sn 104_Elut_KUR_7 \
-xl_sn 77_Elut_KUR_14 \
-xl_sn 143_Elut_KUR_5original \
-xl_sn 144_Elut_KUR_1 \
-xl_sn 145_Elut_KUR_18 \
-xl_sn 154_Elut_KUR_11 \
-xl_sn 155_Elut_KUR_9 \
-xl_sn 90_Elut_KUR_20 \
-xl_sn 80_Elut_KUR_19 \
-xl_sn 91_Elut_KUR_21 \
-xl_sn 92_Elut_KUR_22 \
-xl_sn 93_Elut_KUR_23


python $scriptdir/$script --vcf ${vcfdir}/step-by-step-vcf-sfs/KUR.allCalled.neutralOnly.$vcfFile --pop KUR --outdir ${vcfdir}/step-by-step-vcf-sfs/ --outPREFIX ${vcfFile%_${infile}}


done
