#! /bin/bash
#$ -cwd
#$ -l h_rt=100:00:00,h_data=16G,highp,arch=intel*
#$ -N vcf1a_filtering
#$ -o /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -e /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -m abe
#$ -M ab08028
######## real script starts here:

## Note: when filtering by genotype, the sites that are no call (./.) will get a PASS in the FT field of the genotype
# because they essentially weren't evaluated. This is a really dumb GATK thing. You can fix it with Clare's bespoke script
# But it's also okay, because we filter by the maximum fraction of NO CALL genotypes, NOT by the fraction of PASS genotypes
# So it actually is okay, because the ./. won't contribute to anything

# 20181127: ** big changes for filtering genotypes called 20181119***
# adding an overall per site DP filter that is very loose, min DP 500, no max DP.
# note that GATK doesn't recommend DP filters, esp. for captures, since you expect pileups of reads. 
# these values are arbitrary; based on plotting results on one scaffold for 20180806 genotypes.
# Genotypes: set min genotype DP to 8 instead of 12, with no max DP for genotypes.
# continuing to filter SnpClusters at the same time as HF. I found that doing them sequentially results
# in slightly fewer snpclusters called (~566 for one chromosome) but I want to get rid of sites 
# that are surrounded by low quality snps (within 10bp ), rather than wait until those sites
# have been filtered away and I have no idea that a 'good' site was actually surrounded by multiple bad sites
## Based on plotting the QD dist for one scaffold, see that most sites are clustered with QD > 20, shifted pretty far to the right
# so instead of QD 2 which is from GATK, I am going to switch to QD < 10 as a filter. 
# modules
source /u/local/Modules/default/init/modules.sh
module load java
module load python/2.7
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar
tabix=/u/home/a/ab08028/klohmueldata/annabel_data/bin/tabix-0.2.6/tabix

#### parameters:
rundate=20181119 # date genotypes were called (vcf_20180806 includes capture 02)
noCallFrac=0.2 # maximum fraction of genotypes that can be "no call" (./.) # note that genotypes can still "PASS" if 


#### file locations
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/vcf_filtering
mkdir -p $wd
indir=$SCRATCH/captures/vcfs/vcf_${rundate}
infile=raw_variants.vcf.gz ### make sure this doesn't have a path as part of its name! just infile names
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta

scriptdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_filtering
noCallScript=$scriptdir/filtering_getNoCallPerInd.py

# location of vcf checking and filtering script
# incompatible scaffolds: repeatMaskCoords=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/repeatMaskingCoordinates/masking_coordinates.bed
# repeat masking is optional: my target captrue was designed away from repeat regions, so not a huge deal 
# if you don't include this; but trying to be extra thorough.
# these coords are downloaded from NCBI FTP site for mustela putorius furo ; if using elut genome, use Annotation repeat coords; 
# ftp://ftp.ncbi.nlm.nih.gov/genomes/Mustela_putorius_furo/ <-- location of mustela files
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Mustela_putorius_furo/masking_coordinates.gz

vcfdir=$wd/${rundate}_filtered # date you called genotypes
mkdir -p $vcfdir
mkdir -p $vcfdir/filteringStats

#################################################################################
############################ Prepare File ####################################
################################################################################# 
# trim alternates 
# update: 20180809: stopped filtering by 80% call rate here; want to get rid of bad individuals before I do that.
echo "step 1: trim alternates"
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-V ${indir}/${infile} \
-trimAlternates \
-o ${vcfdir}/'all_1_TrimAlt_'${infile} \
-select "DP > 500"
# made this very lenient, just want to get rid of super crappy sites avg ~8 reads/sample
# below will do further filtering at GT level



#################################################################################
############################ BIALLELIC SNPs ####################################
#################################################################################
# Select only variant sites: (note numbering scheme: 2snp indicates it's step 2 of the snps)
echo "snp step 2: select biallelic snps"
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-V ${vcfdir}/'all_1_TrimAlt_'${infile} \
--restrictAllelesTo BIALLELIC \
--selectTypeToInclude SNP \
-o ${vcfdir}/'snp_2_Filter_TrimAlt_'${infile}


## this:
# 1. does NOT mask out the UCSCRepeats (FAIL_RepMask) (do this later when making SFS)
# 2. applies gatk hard filters (FAIL_GATKHF)
# use genotypeFilterExpression to filter individual genotypes  and **** --setFilteredGtToNocall to change filtered genotype to "no call" (./.) ****
# 3. genotype quality (<20 filtered out) (FAIL_GQ)
# X. *no longer doing this* Individual Depth > 1000 filtered out (FAIL_DP_HIGH)
# 5. Individual DP < 8 filtered out (FAIL_DP_LOW)
# 6. clustered snps (3/10) (SnpCluster) *note, this will filter more snps out than if you did it sequentially with HFs - since it takes the filtered snps into account*
# adding this: --missingValuesInExpressionsShouldEvaluateAsFailing : see how it impacts things
echo "snp step 3: variant filtering"
# modified QD filter to be QD < 10.0 instead of 2.0 based on my plot which showed most variants clustered at QD > 20 (I'm guessing this is a capture feature)
# I separated QD out from the rest of the GATK HFs so that I can better see how it is behaving.
java -jar -Xmx4G ${GATK} \
-T VariantFiltration \
-R ${REFERENCE} \
-V ${vcfdir}/'snp_2_Filter_TrimAlt_'${infile} \
--filterExpression "FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" \
--filterName "FAIL_GATKHF" \
--filterExpression "QD < 10.0" \
--filterName "FAIL_QD10" \
--genotypeFilterExpression "GQ < 20" \
--genotypeFilterName "FAIL_GQ" \
--genotypeFilterExpression "DP < 8" \
--genotypeFilterName "FAIL_DP_LOW" \
--clusterWindowSize 10 --clusterSize 3 \
--setFilteredGtToNocall \
-o ${vcfdir}/'snp_3_Flagged_GQ_DP_GaTKHF_cluster_'${infile}

# removed:
# --genotypeFilterExpression "DP > 1000" \
# --genotypeFilterName "FAIL_DP_HIGH" \

# do variant eval here * **** # 

## note: don't use the 'if it is missing, the site fails' flag: many sites don't have MQRankSum annotations but don't want those to fail
# --mask ${repeatMaskCoords} --maskName "FAIL_RepMask" \
# skipping repeat masking because I designed exome capture away from repeats
# note that if you want to do it, you'll need to make scaffold IDs consistent between the repeat mask (NCBI) and the ref genome
# there's a 1:1 map on the FTP scaffold name document. 

####### PAY ATTENTION TO WARNINGS! IF IT WARNS THAT AN ANNOTATION IS MISSING, YOU'LL HAVE WAY TOO MANY FAILING THAT SHOULDN'T ##########

### Select only passing variants: (this only selects based on those that don't fail the filterExpressions, not the genotypeFilterExpressions)
# 
echo "snp step 4: select passing variants"

java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-V ${vcfdir}/'snp_3_Flagged_GQ_DP_GaTKHF_cluster_'${infile} \
--excludeFiltered \
-trimAlternates \
-o ${vcfdir}/'snp_4_Filtered_GQ_DP_GaTKHF_cluster_'${infile}

#################################################################################
############################ INVARIANT SITES ####################################
#################################################################################
echo "starting: nv step 2: select non variant sites"
# need to redo invariants
## Select the invariants:
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-V ${vcfdir}/'all_1_TrimAlt_'${infile} \
--selectTypeToInclude NO_VARIATION \
--selectTypeToExclude INDEL \
-o ${vcfdir}/'nv_2_AllNonVariants_'${infile}
echo "done: nv step 2: select non variant sites"


echo "starting nv step 3: filter non variant sites"

java -jar -Xmx4G ${GATK} \
-T VariantFiltration \
-R ${REFERENCE} \
-V ${vcfdir}/'nv_2_AllNonVariants_'${infile} \
--filterExpression "QUAL < 30 " \
--filterName "FAIL_QUAL30" \
--genotypeFilterExpression "RGQ < 1" \
--genotypeFilterName "FAIL_RGQ" \
--genotypeFilterExpression "DP < 8" \
--genotypeFilterName "FAIL_DP_LOW" \
--setFilteredGtToNocall \
-o ${vcfdir}/'nv_3_Flagged_DP_RGQ_QUAL_'${infile}

# removing: 
#--genotypeFilterExpression "DP > 1000" \
#--genotypeFilterName "FAIL_DP_HIGH" \

# 20180731 : found bug in my script : was missing --setFilteredGtToNocall for invariant sites. Need to rerun that. 
# took out: --mask ${repeatMaskCoords} --maskName "FAIL_RepMask" \
# note the rgq filter may be redundant: https://software.broadinstitute.org/gatk/blog?id=6495
echo "done nv step 3: filter non variant sites"

### Second round of select variants

echo "starting nv step 4: select only passing non variant sites"

java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-V ${vcfdir}/'nv_3_Flagged_DP_RGQ_QUAL_'${infile} \
--excludeFiltered \
-trimAlternates \
-o ${vcfdir}/'nv_4_Filtered_DP_RGQ_QUAL_'${infile}

# for now taking out: --maxNOCALLfraction $noCallFrac \
## merge back together the variant and invariant files - hopefully we can do that.
#20180731: tyring to add -trimAlternates again at this stage. not sure why it keeps getting through 

echo "done nv step 4: select only passing non variant sites"

#################################################################################
############################ MERGE ALL PASSING SITES ###########################
#################################################################################
echo "starting step 5a: merge passing sites"
java -jar -Xmx4G ${GATK} \
-T CombineVariants \
-R ${REFERENCE} \
-V ${vcfdir}/'snp_4_Filtered_GQ_DP_GaTKHF_cluster_'${infile} \
-V ${vcfdir}/'nv_4_Filtered_DP_RGQ_QUAL_'${infile} \
--assumeIdenticalSamples \
-o ${vcfdir}/'all_5_passingFilters_'${infile}
echo "done step 5a: merge passing sites"


################## Variant Evaluation ############## (still experimental)
# only evaluate on one scaffold (takes too long otherwise)
echo "starting step 5b: variant evaluation"

scaffold="GL896899.1"
java -jar -Xmx4G ${GATK} \
-T VariantEval \
-R $REFERENCE \
-o ${vcfdir}/filteringStats/${scaffold}.'allSNPstages'.variant.eval.txt \
--eval:all1 ${vcfdir}/'all_1_TrimAlt_'${infile} \
--eval:snp2 ${vcfdir}/'snp_2_Filter_TrimAlt_'${infile} \
--eval:snp3 ${vcfdir}/'snp_3_Flagged_GQ_DP_GaTKHF_cluster_'${infile} \
--eval:snp4 ${vcfdir}/'snp_4_Filtered_GQ_DP_GaTKHF_cluster_'${infile} \
-L $scaffold
   


########################## At this stage, calculate the no-call per individual (~10hours) #######################
# takes ~2hrs
echo "starting step 5c: calculate missing calls per individual"
python $noCallScript ${vcfdir}/'all_5_passingFilters_'${infile} ${vcfdir}/filteringStats/'noCall_per_Ind_all_5_passingFilters_'${infile%.vcf.gz}.txt
# then want to remove bad individuals and then run bespoke filters as a final check of AN/AC etc.

echo "done step 5b: calculated missing calls per individual"

