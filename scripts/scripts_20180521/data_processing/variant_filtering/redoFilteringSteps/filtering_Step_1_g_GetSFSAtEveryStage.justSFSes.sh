#! /bin/bash
#$ -cwd
#$ -l h_rt=50:00:00,h_data=16G,highp
#$ -N makeSFS_eachStage
#$ -o /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -e /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -m abe
#$ -M ab08028

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
script=generate1DSFS.py

## bad individuals (make sure gone from early steps?)
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

##############################################

mkdir -p ${vcfdir}/step-by-step-vcf-sfs
infile=raw_variants.vcf.gz
# eventually
for vcfFile in 'all_1_TrimAlt_'${infile} 'all_5_passingFilters_'${infile}
do
# get stats?
echo "starting on " $vcfFile

# Make Commanders SFS
#echo "COM"
#java -jar $GATK \
#-R $REFERENCE \
#-T SelectVariants \
#--variant $vcfdir/$vcfFile \
#-o ${vcfdir}/step-by-step-vcf-sfs/COM.allCalled.neutralOnly.$vcfFile \
#-se '.+_Elut_BER_.+' \
#-se '.+_Elut_MED_.+' \
#--maxNOCALLfraction $perPopNoCallFrac \
#-L $neutralBed

python $scriptdir/$script --vcf ${vcfdir}/step-by-step-vcf-sfs/COM.allCalled.neutralOnly.$vcfFile --pop COM --outdir ${vcfdir}/step-by-step-vcf-sfs/ --outPREFIX ${vcfFile%_${infile}}

# California --> CA (include the RWAB hiseq4000 samples)
#echo "CA"
#java -jar $GATK \
#-R $REFERENCE \
#-T SelectVariants \
#--variant $vcfdir/$vcfFile \
#-o ${vcfdir}/step-by-step-vcf-sfs/CA.allCalled.neutralOnly.$vcfFile \
#-se '.+_Elut_CA_.+' \
#-se 'RWAB003_.+_ELUT_CA_.+' \
#--maxNOCALLfraction $perPopNoCallFrac \
#-L $neutralBed

python $scriptdir/$script --vcf ${vcfdir}/step-by-step-vcf-sfs/CA.allCalled.neutralOnly.$vcfFile --pop CA --outdir ${vcfdir}/step-by-step-vcf-sfs/ --outPREFIX ${vcfFile%_${infile}}


# Alaska --> AK : remove admixed and make pop specific vcf
echo "AK"
#java -jar $GATK \
#-R $REFERENCE \
#-T SelectVariants \
#--variant $vcfdir/$vcfFile \
#-o ${vcfdir}/step-by-step-vcf-sfs/AK.allCalled.neutralOnly.$vcfFile \
#-se '.+_Elut_AK_.+' \
#--maxNOCALLfraction $perPopNoCallFrac \
#-xl_sn ${ind12} \
#-xl_sn ${ind13} \
#-xl_sn ${ind14} \
#-xl_sn ${ind15} \
#-xl_sn ${ind16} \
#-xl_sn ${ind17} \
#-xl_sn ${ind18} \
#-L $neutralBed

python $scriptdir/$script --vcf ${vcfdir}/step-by-step-vcf-sfs/AK.allCalled.neutralOnly.$vcfFile --pop AK --outdir ${vcfdir}/step-by-step-vcf-sfs/ --outPREFIX ${vcfFile%_${infile}}

# Aleutian --> AL
#echo "AL"
#java -jar $GATK \
#-R $REFERENCE \
#-T SelectVariants \
#--variant $vcfdir/$vcfFile \
#-o ${vcfdir}/step-by-step-vcf-sfs/AL.allCalled.neutralOnly.$vcfFile \
#-se '.+_Elut_AL_.+' \
#--maxNOCALLfraction $perPopNoCallFrac \
#-L $neutralBed

python $scriptdir/$script --vcf ${vcfdir}/step-by-step-vcf-sfs/AL.allCalled.neutralOnly.$vcfFile --pop AL --outdir ${vcfdir}/step-by-step-vcf-sfs/ --outPREFIX ${vcfFile%_${infile}}

# Kuril --> KUR (include the RWAB hiseq4000 samples)
echo "KUR"
#java -jar $GATK \
#-R $REFERENCE \
#-T SelectVariants \
#--variant $vcfdir/$vcfFile \
#-o ${vcfdir}/step-by-step-vcf-sfs/KUR.allCalled.neutralOnly.$vcfFile \
#-se '.+_Elut_KUR_.+' \
#-se 'RWAB003_.+_ELUT_KUR_.+' \
#--maxNOCALLfraction $perPopNoCallFrac \
#-xl_sn ${ind7} \
#-xl_sn ${ind8} \
#-xl_sn ${ind9} \
#-xl_sn ${ind10} \
#-xl_sn ${ind11} \
#-L $neutralBed


python $scriptdir/$script --vcf ${vcfdir}/step-by-step-vcf-sfs/KUR.allCalled.neutralOnly.$vcfFile --pop KUR --outdir ${vcfdir}/step-by-step-vcf-sfs/ --outPREFIX ${vcfFile%_${infile}}


done
