#! /bin/bash
#$ -cwd
#$ -l h_rt=5:00:00,h_data=16G
#$ -N testingGATK
#$ -o /u/flashscratch/a/ab08028/sandbox/reports
#$ -e /u/flashscratch/a/ab08028/sandbox/reports
#$ -m abe
#$ -M ab08028

########## testing filtering with elut reference for now
############ test parameters: 1 chromosome; capture 02 data. #############
module load java
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar

numChr=52 # 2 * number of diploids
noCallFrac=0.2 # maximum fraction of genotypes that can be "no call" (./.) # note that genotypes can still "PASS" if 

indir=$SCRATCH/sandbox
wd=$SCRATCH/sandbox

infile=elut.raw_variants.20170914.vcf.gz ### make sure this doesn't have a path as part of its name! just infile name
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/sea_otter_genome/dedup_99_indexed_USETHIS/sea_otter_23May2016_bS9RH.deduped.99.fasta
mychrom=ScbS9RH_6353 # for now  select 1 scaffold for tests, then remove
#repeatMaskCoords=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/repeatMaskingCoordinates/masking_coordinates
outdir=$wd/vcf_filtering
mkdir -p $outdir

###################################################################################################
###################################### run through  1 #############################################
###################################################################################################
# trim alternates and set maxNOCALL to 20% # good
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-L ${mychrom} \
-V ${indir}/${infile} \
-trimAlternates \
--maxNOCALLnumber $maxNoCall \
-o ${outdir}/1_TrimAltRemoveNoCall_${infile}

# testing notes: vcf file requires an index; make sure infile doesn't have path in it.

## biallelic snps:
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-L ${mychrom} \
-V ${outdir}/1_TrimAltRemoveNoCall_${infile} \
--restrictAllelesTo BIALLELIC \
--selectTypeToInclude SNP \
-o ${outdir}/'2snp_Filter_TrimAltRemoveNoCall_'${infile}

# testing notes: no errors



## this applies gatk hard filters, masks out the UCSCRepeats, 
# masks clustered snps (3/10)
java -jar -Xmx4G ${GATK} \
-T VariantFiltration \
-R ${REFERENCE} \
-L ${mychrom} \
-V ${outdir}/'2snp_Filter_TrimAltRemoveNoCall_'${infile} \
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
--missingValuesInExpressionsShouldEvaluateAsFailing \
-o ${outdir}/'3snp_Flagged_VF_GaTKHF_cluster_missingFailing'${infile}
#testing notes: have to remove mask for now (put back in for mfu); no errors


### Annabel's addition: (maybe combine with above)
# want to filter on geotype quality (GQ) and individual site DP (lt 12 for now, gt 1000+?). Maybe up it to 15?
# filtered should be set as no call (./.)
# java -jar -Xmx4G ${GATK} \
# -T VariantFiltration \
# -R ${REFERENCE} \
# -L ${mychrom} \
# -V ${outdir}/'3snp_Flagged_VF_GaTKHF_cluster_'${infile} \
# --genotypeFilterExpression "GQ < 20" \
# --genotypeFilterName "FAIL_gtGQFilter" \
# --setFilteredGtToNocall \
# -o ${outdir}/'4snp_Flagged_GT_VF_GaTKHF_cluster_'${infile}
# testing notes: no errors
################# now at this step: look at what passed and failed
### think about how genotypes have been changed to ./. and if they should have been.
### how to check? 
zcat '2snp_Filter_TrimAltRemoveNoCall_'${infile} | grep -v -c "#" # how many biallelic snps are there 
grep FAIL_vgtFilter '4snp_Flagged_GT_VF_GaTKHF_cluster_'${infile} | less -S 
# seems to be working~ failed individual genotypes and switched them to ./.

###################################################################################################
###################################### run through  2 #############################################
###################################################################################################
# trim alternates and set maxNOCALL to 20% *and do depth filtering per genotype
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-L ${mychrom} \
-V ${outdir}/1_TrimAltRemoveNoCall_${infile} \
--restrictAllelesTo BIALLELIC \
--selectTypeToInclude SNP \
-o ${outdir}/'2snp_Filter_TrimAltRemoveNoCall_'${infile}


## this:
# 1. masks out the UCSCRepeats (FAIL_RepMask)
# 2. applies gatk hard filters (FAIL_GATKHF)
# use genotypeFilterExpression to filter individual genotypes  and **** --setFilteredGtToNocall to change filtered genotype to "no call" (./.) ****
# 3. genotype quality (<20 filtered out) (FAIL_GQ)
# 4. Individual Depth > 1000 filtered out (FAIL_DP_HIGH)
# 5. Individual DP < 12 filtered out (FAIL_DP_LOW)
# 6. clustered snps (3/10) (SnpCluster)

java -jar -Xmx4G ${GATK} \
-T VariantFiltration \
-R ${REFERENCE} \
-L ${mychrom} \
-V ${outdir}/'2snp_Filter_TrimAltRemoveNoCall_'${infile} \
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
-o ${outdir}/'3snp_Flagged_GQ_DP_GaTKHF_cluster_RepMask_'${infile}

#### DOES THIS WORK???? Check very carefully. 
# Does it appropriately update AN? Yes it does.

### Select only passing variants: sites that passed site filters and 
# that have no more than 20% no call sites (no call sites: ./. may be due to failing above genotype filters,
# or may be due to not being called in the first place)
# Note that there is a GATK oddity where some ./. that weren't called in the first place are marked as PASS
# not as fail because they didn't get evaluated. So that is why I'm not filtering based on genotype FT pass/ fail
# instead just by fraction of no call

java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-L ${mychrom} \
-V ${outdir}/'3snp_Flagged_GQ_DP_GaTKHF_cluster_RepMask_'${infile} \
--excludeFiltered \
--maxNOCALLfraction $noCallFrac \
-o ${outdir}/'4snp_Filtered_80percCall_GQ_DP_GaTKHF_cluster_RepMask_'${infile}
# testing notes: appears to work as expected. 

#################################################################################
############################ INVARIANT SITES ####################################
#################################################################################

## Select the invariants - does this work, or should I use jexl to select 'ALT = '.'
# i use the maxnocall fraction because some sites are missing across the board
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-L ${mychrom} \
-V ${outdir}/1_TrimAltRemoveNoCall_${infile} \
--selectTypeToInclude NO_VARIATION \
-o ${outdir}/'2nv_AllNonVariants_'${infile}

# tesing notes: no errors 

java -jar -Xmx4G ${GATK} \
-T VariantFiltration \
-R ${REFERENCE} \
-L ${mychrom} \
-V ${outdir}/'2nv_AllNonVariants_'${infile} \
--filterExpression "QUAL < 30 " \
--filterName "FAIL_QUAL30" \
--genotypeFilterExpression "RGQ < 1" \
--genotypeFilterName "FAIL_RQG" \
--genotypeFilterExpression "DP > 1000" \
--genotypeFilterName "FAIL_DP_HIGH" \
--genotypeFilterExpression "DP < 12" \
--genotypeFilterName "FAIL_DP_LOW" \
-o ${outdir}/'3nv_Flagged_DP_RGQ_QUAL_RepMask_'${infile}
# temporarily have to take out repmask because this isn't mapped to ferret 
# test notes: slow but no errors

### Second round of select variants
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-L ${mychrom} \
-V ${outdir}/'3nv_Flagged_DP_RGQ_QUAL_RepMask_'${infile} \
--excludeFiltered \
--maxNOCALLfraction $noCallFrac \
-o ${outdir}/'4nv_Filtered_80percCall_DP_RGQ_QUAL_RepMask_'${infile}
## merge back together the variant and invariant files - hopefully we can do that.
# no errors;
java -jar -Xmx4G ${GATK} \
-T CombineVariants \
-R ${REFERENCE} \
-V ${outdir}/'4snp_Filtered_80percCall_GQ_DP_GaTKHF_cluster_RepMask_'${infile} \
-V ${outdir}/'4nv_Filtered_80percCall_DP_RGQ_QUAL_RepMask_'${infile} \
--assumeIdenticalSamples \
-o ${outdir}/'5_allSites_passingFilters_80percCall_'${infile}

# select the biallelic snps: 
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-V ${outdir}/'5_allSites_passingFilters_80percCall_'${infile} \
--restrictAllelesTo BIALLELIC \
--selectTypeToInclude SNP \
-o ${outdir}/'snp_5_passingAllFilters_postMerge_'${infile}

## Select the invariants:
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-V ${outdir}/'5_allSites_passingFilters_80percCall_'${infile} \
--selectTypeToInclude NO_VARIATION \
-o ${outdir}/'nv_5_passingAllFilters_postMerge_'${infile}

################## Run through 3: want to try FAILing those that dont' have an entry ####################
