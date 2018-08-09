#! /bin/bash
#$ -cwd
#$ -l h_rt=24:00:00,h_data=16G
#$ -N vcf_filtering
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
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar

#### parameters:
rundate=20180724 # date genotypes were called
noCallFrac=0.2 # maximum fraction of genotypes that can be "no call" (./.) # note that genotypes can still "PASS" if 


#### file locations
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/vcf_filtering
mkdir -p $wd
indir=$SCRATCH/captures/vcfs/vcf_${rundate}
infile=raw_variants.vcf.gz ### make sure this doesn't have a path as part of its name! just infile names
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar
bespokeFilterScript= # location of vcf checking and filtering script
# incompatible scaffolds: repeatMaskCoords=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/repeatMaskingCoordinates/masking_coordinates.bed
# repeat masking is optional: my target captrue was designed away from repeat regions, so not a huge deal 
# if you don't include this; but trying to be extra thorough.
# these coords are downloaded from NCBI FTP site for mustela putorius furo ; if using elut genome, use Annotation repeat coords; 
# ftp://ftp.ncbi.nlm.nih.gov/genomes/Mustela_putorius_furo/ <-- location of mustela files
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Mustela_putorius_furo/masking_coordinates.gz

outdir=$wd/${rundate}_filtered # date you called genotypes
mkdir -p $outdir

#################################################################################
############################ Prepare File ####################################
################################################################################# 
# trim alternates 
# and set maxNOCALL to 20% (have to do this again at the end, but doing it now to make file smaller)
echo "step 1: trim alternates and max no call fraction 20%"
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-V ${indir}/${infile} \
-trimAlternates \
--maxNOCALLfraction $noCallFrac \
-o ${outdir}/'all_1_TrimAlt80Perc_'${infile}



#################################################################################
############################ BIALLELIC SNPs ####################################
#################################################################################
# Select only variant sites: (note numbering scheme: 2snp indicates it's step 2 of the snps)
echo "snp step 2: select biallelic snps"
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-V ${outdir}/'all_1_TrimAlt80Perc_'${infile} \
--restrictAllelesTo BIALLELIC \
--selectTypeToInclude SNP \
-o ${outdir}/'snp_2_Filter_TrimAlt80Perc_'${infile}


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
-V ${outdir}/'snp_2_Filter_TrimAlt80Perc_'${infile} \
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
--maxNOCALLfraction $noCallFrac \
-trimAlternates \
-o ${outdir}/'snp_4_Filtered_80percCall_GQ_DP_GaTKHF_cluster_'${infile}
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
-V ${outdir}/'all_1_TrimAlt80Perc_'${infile} \
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
--maxNOCALLfraction $noCallFrac \
-trimAlternates \
-o ${outdir}/'nv_4_Filtered_80percCall_DP_RGQ_QUAL_'${infile}
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
-V ${outdir}/'snp_4_Filtered_80percCall_GQ_DP_GaTKHF_cluster_'${infile} \
-V ${outdir}/'nv_4_Filtered_80percCall_DP_RGQ_QUAL_'${infile} \
--assumeIdenticalSamples \
-o ${outdir}/'all_5_passingFilters_80percCall_'${infile}
echo "done step 5a: merge passing sites"

#################################################################################
############################ RUN BESPOKE FILTERS and UPDATE AN/AC ##########################
#################################################################################
# These filters will check for:
# 1. ref or alt alleles that aren't a single letter (AGCT) or . (alt)
# 2. genotypes that aren't in 0/0, 0/1, 1/1 or ./. (maybe it's phased, etc)
# 3. must have qual score
# 4. must be PASS for site
# 5. make sure not missing DP, AN, GT, AD, DP or GQ/RGQ
# 6. make sure no called genotype is missing any info from the genotype info field
# 7. gets rid of sites where all calls are 0/1 (all hets)
# 8. updates AN and AC based on final sets of calls (these aren't updated when GATK does genotype filtering)

# this script does NOT: change any genotypes; do any genotype filtering; change any FT fields for genotypes (./. gts will still be PASS if they started as ./. -- bit of GATK weirdness that isn't fatal)
echo "starting step 6a: carrying out bespoke filtering"

python $bespokeFilterScript ${outdir}/'all_5_passingFilters_80percCall_'${infile} ${outdir}/'all_6_passingBespoke_passingFilters_80percCall_'${infile%.gz} ${outdir}/'fail_all_6_FAILINGBespoke_passingFilters_80percCall_'${infile%.vcf.gz}.txt
# gzip the result:
gzip  ${outdir}/'all_6_passingBespoke_passingFilters_80percCall_'${infile%.gz}

echo "done with step 6a: carrying out bespoke filtering"

#################################################################################
############################ Get final sets of variant / invariant###############
#################################################################################
echo "starting step 6b: select final passing snps from merged file"
# Some sites that may have started as variant may have become INVARIANT by the end.
# Want these to end up in the INVARIANT category. 
# select the biallelic snps: 
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-V ${outdir}/'all_6_passingBespoke_passingFilters_80percCall_'${infile} \
--restrictAllelesTo BIALLELIC \
--selectTypeToInclude SNP \
-o ${outdir}/'snp_6_passingBespoke_passingAllFilters_postMerge_'${infile}
echo "done step 6b: select final passing snps from merged file"

## Select the invariants:
echo "starting step 6c: select final passing nonvariant sites from merged file"

java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-V ${outdir}/'all_6_passingBespoke_passingFilters_80percCall_'${infile} \
--selectTypeToInclude NO_VARIATION \
--selectTypeToExclude INDEL \
-o ${outdir}/'nv_6_passingBespoke_passingAllFilters_postMerge_'${infile}
echo "done step 6c: select final passing nonvariant sites from merged file"

##################################################################
############################ GET STATS ###########################
##################################################################
echo "starting getting statistics"
# starting sites:
echo "# Filtering Statistics for SNPs and invariant sites" > ${outdir}/filteringStats.${rundate}.txt
stat0=`zcat ${indir}/${infile} | grep -v -c "#"` 
echo "stat0 starting_sites" $stat0 >> ${outdir}/filteringStats.${rundate}.txt
# sites after trim Altnernates and removing sites with >20% no call
stat1=`zcat ${outdir}/all_1_TrimAlt80Perc_${infile} | grep -v -c "#"`
echo "stat1 sites_afterTrimAlt" $stat1 >> ${outdir}/filteringStats.${rundate}.txt

########## biallelic snps: ##########
stat2=`zcat ${outdir}/'snp_2_Filter_TrimAlt80Perc_'${infile} | grep -v -c "#"`
echo "stat2 biSNPS_initial" $stat2 >> ${outdir}/filteringStats.${rundate}.txt

# snps failing filters:
# note these may overlap (same snp can fail multiple things)

stat31=`zcat ${outdir}/'snp_3_Flagged_GQ_DP_GaTKHF_cluster_'${infile} | grep -v "#" | grep -c FAIL_RepMask`
echo "stat3.1 SNPS_failingREPMASK" $stat31 >> ${outdir}/filteringStats.${rundate}.txt
stat32=`zcat ${outdir}/'snp_3_Flagged_GQ_DP_GaTKHF_cluster_'${infile} | grep -v "#" | grep -c FAIL_GATKHF`
echo "stat3.2 SNPS_failingGATKHF" $stat32 >> ${outdir}/filteringStats.${rundate}.txt
stat36=`zcat ${outdir}/'snp_3_Flagged_GQ_DP_GaTKHF_cluster_'${infile} | grep -v "#" | grep -c SnpCluster`
echo "stat3.6 SNPS_failingSNPCluster" $stat36 >> ${outdir}/filteringStats.${rundate}.txt
# genotypes that fail:
# note that grep -o outputs ALL occurences, not all lines (so you can count multiple on one line)
stat33=`zcat ${outdir}/'snp_3_Flagged_GQ_DP_GaTKHF_cluster_'${infile} | grep -v "#" | grep -o FAIL_GQ | wc -l`
echo "stat3.3 var_GTs_failingGQ" $stat33 >> ${outdir}/filteringStats.${rundate}.txt
stat34=`zcat ${outdir}/'snp_3_Flagged_GQ_DP_GaTKHF_cluster_'${infile} | grep -v "#" | grep -o FAIL_DP_HIGH | wc -l`
echo "stat3.4 var_GTs_failingDP_HIGH" $stat34 >> ${outdir}/filteringStats.${rundate}.txt
stat35=`zcat ${outdir}/'snp_3_Flagged_GQ_DP_GaTKHF_cluster_'${infile} | grep -v "#" | grep -o FAIL_DP_LOW| wc -l`
echo "stat3.5 var_GTs_failingDP_LOW" $stat35 >> ${outdir}/filteringStats.${rundate}.txt

# total passing snps
stat4=`zcat ${outdir}/'snp_4_Filtered_80percCall_GQ_DP_GaTKHF_cluster_'${infile} | grep -v -c "#"`
echo "stat4 totalPassingSNPsSites_80Perc_preMerge" $stat4 >> ${outdir}/filteringStats.${rundate}.txt

########## invariant sites: ###########
stat5=`zcat ${outdir}/'nv_2_AllNonVariants_'${infile} | grep -v -c "#"` 
echo "stat5 invarSites_initial" $stat5 >> ${outdir}/filteringStats.${rundate}.txt

# filters of invariant sites:
stat61=`zcat ${outdir}/'nv_3_Flagged_DP_RGQ_QUAL_'${infile} | grep -v "#" | grep -c FAIL_RepMask`
echo "stat6.1 invarSites_failingREPMASK" $stat61 >> ${outdir}/filteringStats.${rundate}.txt
stat62=`zcat ${outdir}/'nv_3_Flagged_DP_RGQ_QUAL_'${infile} | grep -v "#" | grep -c FAIL_QUAL30`
echo "stat6.2 invarSites_failingQUAL30" $stat62 >> ${outdir}/filteringStats.${rundate}.txt
# invariant genotype filters:
stat63=`zcat ${outdir}/'nv_3_Flagged_DP_RGQ_QUAL_'${infile} | grep -v "#" | grep -o FAIL_DP_LOW | wc -l`
echo "stat6.3 invar_GTs_failingDP_LOW" $stat63 >> ${outdir}/filteringStats.${rundate}.txt
stat64=`zcat ${outdir}/'nv_3_Flagged_DP_RGQ_QUAL_'${infile} | grep -v "#" | grep -o FAIL_DP_HIGH | wc -l`
echo "stat6.4 invar_GTs_failingDP_HIGH" $stat64 >> ${outdir}/filteringStats.${rundate}.txt
stat65=`zcat ${outdir}/'nv_3_Flagged_DP_RGQ_QUAL_'${infile} | grep -v "#" | grep -o FAIL_RGQ | wc -l`
echo "stat6.5 invar_GTs_failingRGQ" $stat65 >> ${outdir}/filteringStats.${rundate}.txt

# total passing invariant sites 
stat7=`zcat ${outdir}/'nv_4_Filtered_80percCall_DP_RGQ_QUAL_'${infile} | grep -v -c "#"`
echo "stat7 total_invarSites_passing_80Perc_preMerge" $stat7 >> ${outdir}/filteringStats.${rundate}.txt

################## merged sites ##############
# total passing sites
stat81=`zcat ${outdir}/'all_5_passingFilters_80percCall_'${infile} | grep -v -c "#"`
echo "stat8.1 totalPassingSites_all_preBespoke" $stat81 >> ${outdir}/filteringStats.${rundate}.txt

stat82=`zcat ${outdir}/'all_6_passingBespoke_passingFilters_80percCall_'${infile} | grep -v -c "#"`
echo "stat8.2 totalPassingSites_all_postBespoke" $stat82 >> ${outdir}/filteringStats.${rundate}.txt

stat83=`grep -v -c "#" ${outdir}/'fail_all_6_FAILINGBespoke_passingFilters_80percCall_'${infile%.vcf.gz}.txt`
echo "stat8.3 sites_failing_bespoke" $stat83 >> ${outdir}/filteringStats.${rundate}.txt

# note that some snps may have become invariant after filtering.
stat9=`zcat ${outdir}/'snp_6_passingBespoke_passingAllFilters_postMerge_'${infile} | grep -v -c "#"` 
echo "stat9 totalPassingSNPsSites_80Perc_postMerge_postBespoke" $stat9 >> ${outdir}/filteringStats.${rundate}.txt

stat10=`zcat ${outdir}/'nv_6_passingBespoke_passingAllFilters_postMerge_'${infile} | grep -v -c "#"` 
echo "stat10 total_invarSites_passing_80Perc_postMerge_postBespoke" $stat10 >> ${outdir}/filteringStats.${rundate}.txt
echo "done getting statistics"

## optional : clean up
#rm ${outdir}/1_TrimAltRemoveNoCall_${infile}
#rm ${outdir}/'2snp_Filter_TrimAltRemoveNoCall_'${infile}
#rm ${outdir}/'3snp_VF_Flagged_GaTKHFGQ20_cluster_'${infile}
#rm ${outdir}/'4snp_VF_Filtered_GaTKHFGQ20_cluster_'${infile}
#rm ${outdir}/2nv_AllNonVariants_${infile}
#rm ${outdir}/'3nv_Flagged_AllNonVariantsRepeatmaskQUAL_'${infile}
#rm ${outdir}/'4nv_Filtered_AllNonVariantsRepeatmaskQUAL_'${infile}


