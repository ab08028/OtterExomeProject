#! /bin/bash
#$ -cwd
#$ -l h_rt=24:00:00,h_data=16G
#$ -N vcf_filtering
#$ -o /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -e /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -m abe
#$ -M ab08028
######## real script starts here:
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

repeatMaskCoords=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/repeatMaskingCoordinates/masking_coordinates
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
java -jar -Xmx4G ${GATK} \
-T VariantFiltration \
-R ${REFERENCE} \
-V ${outdir}/'snp_2_Filter_TrimAlt80Perc_'${infile} \
--mask ${repeatMaskCoords} --maskName "FAIL_RepMask" \
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
-o ${outdir}/'snp_3_Flagged_GQ_DP_GaTKHF_cluster_RepMask_'${infile}

####### PAY ATTENTION TO WARNINGS! IF IT WARNS THAT AN ANNOTATION IS MISSING, YOU'LL HAVE WAY TOO MANY FAILING THAT SHOULDN'T ##########

### Select only passing variants: (this only selects based on those that don't fail the filterExpressions, not the genotypeFilterExpressions)
# and only select sites where max 20% of genotypes are not called
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-V ${outdir}/'snp_3_Flagged_GQ_DP_GaTKHF_cluster_RepMask_'${infile} \
--excludeFiltered \
--maxNOCALLfraction $noCallFrac \
-o ${outdir}/'snp_4_Filtered_80percCall_GQ_DP_GaTKHF_cluster_RepMask_'${infile}

### Also want to restrict for sites with calls in 80% of chromosomes. (use AN =?)
### Also want to eliminate any that are all hets. How to find those?
#################################################################################
############################ INVARIANT SITES ####################################
#################################################################################

## Select the invariants:
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-V ${outdir}/1_TrimAltRemoveNoCall_${infile} \
--selectTypeToInclude NO_VARIATION \
-o ${outdir}/'nv_2_AllNonVariants_'${infile}

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
--mask ${repeatMaskCoords} --maskName "FAIL_RepMask" \
-o ${outdir}/'nv_3_Flagged_DP_RGQ_QUAL_RepMask_'${infile}
# 
### Second round of select variants
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-V ${outdir}/'nv_3_Flagged_DP_RGQ_QUAL_RepMask_'${infile} \
--excludeFiltered \
--maxNOCALLfraction $noCallFrac \
-o ${outdir}/'nv_4_Filtered_80percCall_DP_RGQ_QUAL_RepMask_'${infile}
## merge back together the variant and invariant files - hopefully we can do that.

#################################################################################
############################ MERGE ALL PASSING SITES ###########################
#################################################################################
java -jar -Xmx4G ${GATK} \
-T CombineVariants \
-R ${REFERENCE} \
-V ${outdir}/'snp_4_Filtered_80percCall_GQ_DP_GaTKHF_cluster_RepMask_'${infile} \
-V ${outdir}/'nv_4_Filtered_80percCall_DP_RGQ_QUAL_RepMask_'${infile} \
--assumeIdenticalSamples \
-o ${outdir}/'all_5_passingFilters_80percCall_'${infile}


#################################################################################
############################ Get final sets of variant / invariant###############
#################################################################################

# Some sites that may have started as variant may have become INVARIANT by the end.
# Want these to end up in the INVARIANT category. 
# select the biallelic snps: 
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-V ${outdir}/'all_5_passingFilters_80percCall_'${infile} \
--restrictAllelesTo BIALLELIC \
--selectTypeToInclude SNP \
-o ${outdir}/'snp_5_passingAllFilters_postMerge_'${infile}

## Select the invariants:
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-V ${outdir}/'all_5_passingFilters_80percCall_'${infile} \
--selectTypeToInclude NO_VARIATION \
-o ${outdir}/'nv_5_passingAllFilters_postMerge_'${infile}

##################################################################
############################ GET STATS ###########################
##################################################################
# starting sites:
echo "# Filtering Statistics for SNPs and invariant sites" > ${outdir}/filteringStats.${rundate}.txt
stat0=`zcat ${indir}/${infile} | grep -v -c "#"` 
echo "stat0 starting_sites" $stat0 >> ${outdir}/filteringStats.${rundate}.txt
# sites after trim Altnernates and removing sites with >20% no call
stat1=`zcat ${outdir}/all_1_TrimAlt80Perc_${infile} | grep -v -c "#"`
echo "stat1 sites_afterTrimAlt" $stat1 >> ${outdir}/filteringStats.${rundate}.txt

########## biallelic snps: ##########
stat2=`zcat 'snp_2_Filter_TrimAlt80Perc_'${infile} | grep -v -c "#"`
echo "stat2 biSNPS_initial" $stat2 >> ${outdir}/filteringStats.${rundate}.txt

# snps failing filters:
# note these may overlap (same snp can fail multiple things)

stat31=`zcat ${outdir}/'snp_3_Flagged_GQ_DP_GaTKHF_cluster_RepMask_'${infile} | grep -v "#" | grep -c FAIL_RepMask`
echo "stat3.1 SNPS_failingREPMASK" $stat31 >> ${outdir}/filteringStats.${rundate}.txt
stat32=`zcat ${outdir}/'snp_3_Flagged_GQ_DP_GaTKHF_cluster_RepMask_'${infile} | grep -v "#" | grep -c FAIL_GATKHF`
echo "stat3.2 SNPS_failingGATKHF" $stat32 >> ${outdir}/filteringStats.${rundate}.txt
stat36=`zcat ${outdir}/'snp_3_Flagged_GQ_DP_GaTKHF_cluster_RepMask_'${infile} | grep -v "#" | grep -c SnpCluster`
echo "stat3.6 SNPS_failingSNPCluster" $stat36 >> ${outdir}/filteringStats.${rundate}.txt
# genotypes that fail:
# note that grep -o outputs ALL occurences, not all lines (so you can count multiple on one line)
stat33=`zcat ${outdir}/'snp_3_Flagged_GQ_DP_GaTKHF_cluster_RepMask_'${infile} | grep -v "#" | grep -o FAIL_GQ | wc -l`
echo "stat3.3 var_GTs_failingGQ" $stat33 >> ${outdir}/filteringStats.${rundate}.txt
stat34=`zcat ${outdir}/'snp_3_Flagged_GQ_DP_GaTKHF_cluster_RepMask_'${infile} | grep -v "#" | grep -o FAIL_DP_HIGH | wc -l`
echo "stat3.4 var_GTs_failingDP_HIGH" $stat34 >> ${outdir}/filteringStats.${rundate}.txt
stat35=`zcat ${outdir}/'snp_3_Flagged_GQ_DP_GaTKHF_cluster_RepMask_'${infile} | grep -v "#" | grep -o FAIL_DP_LOW| wc -l`
echo "stat3.5 var_GTs_failingDP_LOW" $stat35 >> ${outdir}/filteringStats.${rundate}.txt

# total passing snps
stat4=`zcat 'snp_4_Filtered_80percCall_GQ_DP_GaTKHF_cluster_RepMask_'${infile} | grep -v -c "#"`
echo "stat4 totalPassingSNPsSites_80Perc_preMerge" $stat4 >> ${outdir}/filteringStats.${rundate}.txt

########## invariant sites: ###########
stat5=`zcat 'nv_2_AllNonVariants_'${infile} | grep -v -c "#"` 
echo "stat5 invarSites_initial" $stat5 >> ${outdir}/filteringStats.${rundate}.txt

# filters of invariant sites:
stat61=`zcat ${outdir}/'nv_3_Flagged_DP_RGQ_QUAL_RepMask_'${infile} | grep -v "#" | grep -c FAIL_RepMask`
echo "stat6.1 invarSites_failingREPMASK" $stat61 >> ${outdir}/filteringStats.${rundate}.txt
stat62=`zcat ${outdir}/'nv_3_Flagged_DP_RGQ_QUAL_RepMask_'${infile} | grep -v "#" | grep -c FAIL_QUAL30`
echo "stat6.2 invarSites_failingQUAL30" $stat62 >> ${outdir}/filteringStats.${rundate}.txt
# invariant genotype filters:
stat63=`zcat ${outdir}/'nv_3_Flagged_DP_RGQ_QUAL_RepMask_'${infile} | grep -v "#" | grep -o FAIL_DP_LOW | wc -l`
echo "stat6.3 invar_GTs_failingDP_LOW" $stat63 >> ${outdir}/filteringStats.${rundate}.txt
stat64=`zcat ${outdir}/'nv_3_Flagged_DP_RGQ_QUAL_RepMask_'${infile} | grep -v "#" | grep -o FAIL_DP_HIGH | wc -l`
echo "stat6.4 invar_GTs_failingDP_HIGH" $stat64 >> ${outdir}/filteringStats.${rundate}.txt
stat65=`zcat ${outdir}/'nv_3_Flagged_DP_RGQ_QUAL_RepMask_'${infile} | grep -v "#" | grep -o FAIL_RGQ | wc -l`
echo "stat6.5 invar_GTs_failingRGQ" $stat65 >> ${outdir}/filteringStats.${rundate}.txt

# total passing invariant sites 
stat7=`zcat 'nv_4_Filtered_80percCall_DP_RGQ_QUAL_RepMask_'${infile} | grep -v -c "#"`
echo "stat7 total_invarSites_passing_80Perc_preMerge" $stat7 >> ${outdir}/filteringStats.${rundate}.txt

################## merged sites ##############
# total passing sites
stat8=`zcat 'all_5_passingFilters_80percCall_'${infile} | grep -v -c "#"`
echo "stat8 totalPassingSites_all" $stat8 >> ${outdir}/filteringStats.${rundate}.txt
# note that some snps may have become invariant after filtering.
stat9=`zcat 'snp_5_passingAllFilters_postMerge_'${infile} | grep -v -c "#"` 
echo "stat9 totalPassingSNPsSites_80Perc_postMerge" $stat9 >> ${outdir}/filteringStats.${rundate}.txt

stat10=`zcat 'nv_5_passingAllFilters_postMerge_'${infile} | grep -v -c "#"` 
echo "stat10 total_invarSites_passing_80Perc_postMerge" $stat10 >> ${outdir}/filteringStats.${rundate}.txt

## optional : clean up
#rm ${outdir}/1_TrimAltRemoveNoCall_${infile}
#rm ${outdir}/'2snp_Filter_TrimAltRemoveNoCall_'${infile}
#rm ${outdir}/'3snp_VF_Flagged_GaTKHFGQ20_cluster_'${infile}
#rm ${outdir}/'4snp_VF_Filtered_GaTKHFGQ20_cluster_'${infile}
#rm ${outdir}/2nv_AllNonVariants_${infile}
#rm ${outdir}/'3nv_Flagged_AllNonVariantsRepeatmaskQUAL_'${infile}
#rm ${outdir}/'4nv_Filtered_AllNonVariantsRepeatmaskQUAL_'${infile}


