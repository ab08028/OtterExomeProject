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

# modules
source /u/local/Modules/default/init/modules.sh
module load java
module load python/2.7
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar
tabix=/u/home/a/ab08028/klohmueldata/annabel_data/bin/tabix-0.2.6/tabix

#### parameters:
rundate=20180806 # date genotypes were called (vcf_20180806 includes capture 02)
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

outdir=$wd/${rundate}_filtered # date you called genotypes
mkdir -p $outdir
mkdir -p $outdir/filteringStats

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
-o ${outdir}/'all_1_TrimAlt_'${infile} \
--maxNOCALLfraction 0.5

# haven't run it this way yet, but am getting rid of sites where 50% of sites are no-call TO MAKE THE FILE SMALLER. 
# [later in the steps: will do final filtering to restrict to sites where max 20% of sites are nocall]
# previously when I did it, there were sites with only 1-2 calls that inflate file sizes and slow stuff down
# I wanted to not filter on missingness yet, because I am going to be removing bad individuals that drag down everybody's missingness %%
# but I think there's a middle ground perhaps. At this stage I could filter sites that have >80% missing. (Whereas later I filter more stringently with those 
# that have >20% missing)


#################################################################################
############################ BIALLELIC SNPs ####################################
#################################################################################
# Select only variant sites: (note numbering scheme: 2snp indicates it's step 2 of the snps)
echo "snp step 2: select biallelic snps"
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-V ${outdir}/'all_1_TrimAlt_'${infile} \
--restrictAllelesTo BIALLELIC \
--selectTypeToInclude SNP \
-o ${outdir}/'snp_2_Filter_TrimAlt_'${infile}


## this:
# 1. masks out the UCSCRepeats (FAIL_RepMask)
# 2. applies gatk hard filters (FAIL_GATKHF)
# use genotypeFilterExpression to filter individual genotypes  and **** --setFilteredGtToNocall to change filtered genotype to "no call" (./.) ****
# 3. genotype quality (<20 filtered out) (FAIL_GQ)
# 4. Individual Depth > 1000 filtered out (FAIL_DP_HIGH)
# 5. Individual DP < 12 filtered out (FAIL_DP_LOW)
# 6. clustered snps (3/10) (SnpCluster)
# adding this: --missingValuesInExpressionsShouldEvaluateAsFailing : see how it impacts things
echo "snp step 3: variant filtering"

java -jar -Xmx4G ${GATK} \
-T VariantFiltration \
-R ${REFERENCE} \
-V ${outdir}/'snp_2_Filter_TrimAlt_'${infile} \
--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" \
--filterName "FAIL_GATKHF" \
--genotypeFilterExpression "GQ < 20" \
--genotypeFilterName "FAIL_GQ" \
--genotypeFilterExpression "DP > 1000" \
--genotypeFilterName "FAIL_DP_HIGH" \
--genotypeFilterExpression "DP < 12" \
--genotypeFilterName "FAIL_DP_LOW" \
--clusterWindowSize 10 --clusterSize 3 \
--setFilteredGtToNocall \
-o ${outdir}/'snp_3_Flagged_GQ_DP_GaTKHF_cluster_'${infile}
## note: don't use the 'if it is missing, the site fails' flag: manny sites don't have MQRankSum annotations but don't want those to fail
# --mask ${repeatMaskCoords} --maskName "FAIL_RepMask" \
# skipping repeat masking because I designed exome capture away from repeats
# note that if you want to do it, you'll need to make scaffold IDs consistent between the repeat mask (NCBI) and the ref genome
# there's a 1:1 map on the FTP scaffold name document. 

####### PAY ATTENTION TO WARNINGS! IF IT WARNS THAT AN ANNOTATION IS MISSING, YOU'LL HAVE WAY TOO MANY FAILING THAT SHOULDN'T ##########

### Select only passing variants: (this only selects based on those that don't fail the filterExpressions, not the genotypeFilterExpressions)
# and only select sites where max 20% of genotypes are not called
echo "snp step 4: select passing variants"

java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-V ${outdir}/'snp_3_Flagged_GQ_DP_GaTKHF_cluster_'${infile} \
--excludeFiltered \
-trimAlternates \
-o ${outdir}/'snp_4_Filtered_GQ_DP_GaTKHF_cluster_'${infile}
# for now removing --maxNOCALLfraction $noCallFrac \
# also adding trimAlternates again in case some got through [luckily none got through]
### Also want to restrict for sites with calls in 80% of chromosomes. (use AN =?)
### Also want to eliminate any that are all hets. How to find those?
#################################################################################
############################ INVARIANT SITES ####################################
#################################################################################
echo "starting: nv step 2: select non variant sites"
# need to redo invariants
## Select the invariants:
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-V ${outdir}/'all_1_TrimAlt_'${infile} \
--selectTypeToInclude NO_VARIATION \
--selectTypeToExclude INDEL \
-o ${outdir}/'nv_2_AllNonVariants_'${infile}
echo "done: nv step 2: select non variant sites"


echo "starting nv step 3: filter non variant sites"

java -jar -Xmx4G ${GATK} \
-T VariantFiltration \
-R ${REFERENCE} \
-V ${outdir}/'nv_2_AllNonVariants_'${infile} \
--filterExpression "QUAL < 30 " \
--filterName "FAIL_QUAL30" \
--genotypeFilterExpression "RGQ < 1" \
--genotypeFilterName "FAIL_RGQ" \
--genotypeFilterExpression "DP > 1000" \
--genotypeFilterName "FAIL_DP_HIGH" \
--genotypeFilterExpression "DP < 12" \
--genotypeFilterName "FAIL_DP_LOW" \
--setFilteredGtToNocall \
-o ${outdir}/'nv_3_Flagged_DP_RGQ_QUAL_'${infile}
# 20180731 : found bug in my script : was missing --setFilteredGtToNocall for invariant sites. Need to rerun that. 
# took out: --mask ${repeatMaskCoords} --maskName "FAIL_RepMask" \
# note the rgq filter may be redundant: https://software.broadinstitute.org/gatk/blog?id=6495
echo "done nv step 3: filter non variant sites"

### Second round of select variants

echo "starting nv step 4: select only passing non variant sites"

java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-V ${outdir}/'nv_3_Flagged_DP_RGQ_QUAL_'${infile} \
--excludeFiltered \
-trimAlternates \
-o ${outdir}/'nv_4_Filtered_DP_RGQ_QUAL_'${infile}
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
-V ${outdir}/'snp_4_Filtered_GQ_DP_GaTKHF_cluster_'${infile} \
-V ${outdir}/'nv_4_Filtered_DP_RGQ_QUAL_'${infile} \
--assumeIdenticalSamples \
-o ${outdir}/'all_5_passingFilters_'${infile}
echo "done step 5a: merge passing sites"


########################## At this stage, calculate the no-call per individual (~10hours) #######################
# takes ~2hrs
echo "starting step 5b: calculate missing calls per individual"
python $noCallScript ${outdir}/'all_5_passingFilters_'${infile} ${outdir}/filteringStats/'noCall_per_Ind_all_5_passingFilters_'${infile%.vcf.gz}.txt
# then want to remove bad individuals and then run bespoke filters as a final check of AN/AC etc.

echo "done step 5b: calculated missing calls per individual"

