
######## real script starts here:
# get 80% value
numInd=122 # diploids
# 80% 
###### CAREFUL -- set these manually for now (automate later) ########
perc80=97 # 122 * .8; round down
maxNoCall=25 # 122 - 97  # round up
rundate= # date genotypes were called


# modules
source /u/local/Modules/default/init/modules.sh
module load java
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar

#### file locations
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/vcf_filtering
mkdir -p $wd
vcfdir=$SCRATCH/captures/vcfs
infile=$vcfdir/raw_variants.vcf.gz
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
# trim alternates and set maxNOCALL to 20% # good
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-L ${mychrom} \
-V ${indir}/${infile} \
-trimAlternates \
--maxNOCALLnumber $maxNoCall \
-o ${outdir}/1_TrimAltRemoveNoCall_${infile}

# testing notes: 

#################################################################################
############################ BIALLELIC SNPs ####################################
#################################################################################
# Select only variant sites: (note numbering scheme: 2snp indicates it's step 2 of the snps)

java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-L ${mychrom} \
-V ${outdir}/1_TrimAltRemoveNoCall_${infile} \
--restrictAllelesTo BIALLELIC \
--selectTypeToInclude SNP \
-o ${outdir}/'2snp_Filter_TrimAltRemoveNoCall_'${infile}


## this applies gatk hard filters, masks out the UCSCRepeats
# masks clustered snps (3/10)

java -jar -Xmx4G ${GATK} \
-T VariantFiltration \
-R ${REFERENCE} \
-L ${mychrom} \
-V ${outdir}/'2snp_Filter_TrimAltRemoveNoCall_'${infile} \
--mask ${repeatMaskCoords} --maskName "FAIL_RepMask" \
--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" \
--filterName "FAIL_GATKHF" \
--clusterWindowSize 10 --clusterSize 3 \
-o ${outdir}/'3snp_Flagged_VF_GaTKHF_cluster_'${infile}

### Annabel's addition: (maybe combine with above)
# want to filter on geotype quality (GQ) and individual site DP (lt 12 for now, gt 1000+?). Maybe up it to 15?
# filtered should be set as no call (./.)
java -jar -Xmx4G ${GATK} \
-T VariantFiltration \
-R ${REFERENCE} \
-L ${mychrom} \
-V ${outdir}/'3snp_Flagged_VF_GaTKHF_cluster_'${infile} \
--genotypeFilterExpression "GQ < 20 || DP < 12 || DP > 1000" \
--genotypeFilterName "FAIL_vgtFilter" \
--setFilteredGtToNocall \
-o ${outdir}/'4snp_Flagged_GT_VF_GaTKHF_cluster_'${infile}
#### DOES THIS WORK???? Check very carefully. 

### Select only passing variants:
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-L ${mychrom} \
-V ${outdir}/'4snp_Flagged_GT_VF_GaTKHF_cluster_'${infile} \
-ef \
-o ${outdir}/'5snp_Filtered_GT_VF_GaTKHF_cluster_'${infile}


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
-o ${outdir}/2nv_AllNonVariants_${infile}

java -jar -Xmx4G ${GATK} \
-T VariantFiltration \
-R ${REFERENCE} \
-L ${mychrom} \
-V ${outdir}/2nv_AllNonVariants_${infile} \
--filterExpression "QUAL < 30 " \
--filterName "FAIL_QUALlt30" \
--mask ${repeatMaskCoords} --maskName "FAIL_RepMask" \
-o ${outdir}/'3nv_Flagged_AllNonVariantsRepeatmaskQUAL_'${infile}

### annabel's addition: filter by genotype GQ and Depth
java -jar -Xmx4G ${GATK} \
-T VariantFiltration \
-R ${REFERENCE} \
-L ${mychrom} \
-V ${outdir}/'3nv_Flagged_AllNonVariantsRepeatmaskQUAL_'${infile} \
--genotypeFilterExpression "RGQ < 1 || DP < 12 || DP > 1000" \
--genotypeFilterName "FAIL_nvgtFilter" \
--setFilteredGtToNocall \
-o ${outdir}/'4nv_Flagged_GT_AllNonVariantsRepeatmaskQUAL_'${infile}

### Second round of select variants
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-L ${mychrom} \
-V ${outdir}/'4nv_Flagged_GT_AllNonVariantsRepeatmaskQUAL_'${infile} \
-ef \
-o ${outdir}/'5nv_Filtered_GT_AllNonVariantsRepeatmaskQUAL_'${infile}
## merge back together the variant and invariant files - hopefully we can do that.

#################################################################################
############################ ALL PASSING SITES MERGED ####################################
#################################################################################
java -jar -Xmx4G ${GATK} \
-T CombineVariants \
-R ${REFERENCE} \
-V ${outdir}/'4snp_VF_Filtered_GaTKHF_cluster_'${infile} \
-V ${outdir}/'5nv_Filtered_GT_AllNonVariantsRepeatmaskQUAL_'${infile} \
-L ${mychrom} \
--assumeIdenticalSamples \
-o ${outdir}/'5_allSites_passingFilters__'${infile}

# get stats before you clean update

#step1ct=`grep -v -c "#" ${outdir}/1_TrimAltRemoveNoCall_${infile}`
#step1ct=`grep -v -c "#" ${outdir}/1_TrimAltRemoveNoCall_${infile}`
#step1ct=`grep -v -c "#" ${outdir}/1_TrimAltRemoveNoCall_${infile}`
#step1ct=`grep -v -c "#" ${outdir}/1_TrimAltRemoveNoCall_${infile}`
#step1ct=`grep -v -c "#" ${outdir}/1_TrimAltRemoveNoCall_${infile}`
#echo "1_TrimAltRemoveNoCall $step1ct" >> filteringStats.txt
#echo "1_TrimAltRemoveNoCall $step1ct" >> filteringStats.txt

## clean up
#rm ${outdir}/1_TrimAltRemoveNoCall_${infile}
#rm ${outdir}/'2snp_Filter_TrimAltRemoveNoCall_'${infile}
#rm ${outdir}/'3snp_VF_Flagged_GaTKHFGQ20_cluster_'${infile}
#rm ${outdir}/'4snp_VF_Filtered_GaTKHFGQ20_cluster_'${infile}
#rm ${outdir}/2nv_AllNonVariants_${infile}
#rm ${outdir}/'3nv_Flagged_AllNonVariantsRepeatmaskQUAL_'${infile}
#rm ${outdir}/'4nv_Filtered_AllNonVariantsRepeatmaskQUAL_'${infile}


