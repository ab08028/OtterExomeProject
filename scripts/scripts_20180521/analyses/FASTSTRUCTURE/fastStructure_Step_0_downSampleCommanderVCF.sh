#! /bin/bash
#$ -cwd
#$ -l h_rt=50:00:00,h_data=16G,highp
#$ -N downSampleCommanders_forStructure
#$ -o /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -e /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -m abe
#$ -M ab08028
# this can be run in the shell?
# modules
source /u/local/Modules/default/init/modules.sh
module load java
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar

#### parameters:
rundate=20180806 # date genotypes were called (vcf_20180806 includes capture 02)
noCallFrac=0.2 # maximum fraction of genotypes that can be "no call" (./.) # note that genotypes can still "PASS" if 


#### file locations
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/vcf_filtering/
mkdir -p $wd
infile=raw_variants.vcf.gz ### make sure this doesn't have a path as part of its name! just infile names
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta
vcfdir=$wd/${rundate}_filtered # date you called genotypes

##### Downsample Commanders for faststructure (randomly exlcuding these 21   ; leaving 20 in)
# excluding:
ind1=130_Elut_BER_97
ind2=133_Elut_BER_37
ind3=135_Elut_BER_42
ind4=137_Elut_BER_88
ind5=156_Elut_MED_27
ind6=158_Elut_MED_31
ind7=162_Elut_MED_21
ind8=47_Elut_BER_32
ind9=49_Elut_BER_36
ind10=51_Elut_BER_39
ind11=53_Elut_BER_93
ind12=59_Elut_MED_13
ind13=61_Elut_MED_15
ind14=63_Elut_MED_17
ind15=71_Elut_BER_44
ind16=73_Elut_BER_50
ind17=82_Elut_MED_20
ind18=84_Elut_MED_23
ind19=85_Elut_MED_24
ind20=87_Elut_MED_26
ind21=88_Elut_MED_29
# randomly exclude 50% of MED and 50% of BER
mkdir -p $vcfdir/downsampledVCFs/

# remove from snp file:
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/'snp_7_maxNoCallFrac_'${noCallFrac}'_passingBespoke_passingAllFilters_postMerge_'${infile} \
-o ${vcfdir}/downsampledVCFs/'snp_7_downSampCOM_maxNoCallFrac_'${noCallFrac}'_passingBespoke_passingAllFilters_postMerge_'${infile} \
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
-xl_sn ${ind15} \
-xl_sn ${ind16} \
-xl_sn ${ind17} \
-xl_sn ${ind18} \
-xl_sn ${ind19} \
-xl_sn ${ind20} \
-xl_sn ${ind21}