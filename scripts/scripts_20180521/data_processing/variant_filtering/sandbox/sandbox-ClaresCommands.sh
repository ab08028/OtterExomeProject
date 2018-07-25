# clare step 1:
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-L ${mychrom} \
-V ${indir}/${infile} \
-trimAlternates \
--maxNOCALLnumber 2 \
-o ${outdir}/1_TrimAltRemoveNoCall_${infile}

# ## Get just the bi SNP

java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-L ${mychrom} \
-V ${outdir}/1_TrimAltRemoveNoCall_${infile} \
--restrictAllelesTo BIALLELIC \
--selectTypeToInclude SNP \
-o ${outdir}/'2snp_Filter_TrimAltRemoveNoCall_'${infile}


## this applies gatk hard filters, masks out the UCSCRepeats, and makes genotypes with GQ <20 or DP <2 > max ... 

java -jar -Xmx4G ${GATK} \
-T VariantFiltration \
-R ${REFERENCE} \
-L ${mychrom} \
-V ${outdir}/'2snp_Filter_TrimAltRemoveNoCall_'${infile} \
--mask ${UCSCRepeats} --maskName "FAIL_RepMask" \
--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" \
--filterName "FAIL_GATKHF" \
--clusterWindowSize 10 --clusterSize 3 \
-o ${outdir}/'3snp_VF_Flagged_GaTKHF_cluster_'${infile}

### Second round of select variants
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-L ${mychrom} \
-V ${outdir}/'3snp_VF_Flagged_GaTKHF_cluster_'${infile} \
-ef \
-o ${outdir}/'4snp_VF_Filtered_GaTKHF_cluster_'${infile}

## Get just the nonvariants

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
--mask ${UCSCRepeats} --maskName "FAIL_RepMask" \
-o ${outdir}/'3nv_Flagged_AllNonVariantsRepeatmaskQUAL_'${infile}

### Second round of select variants
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-L ${mychrom} \
-V ${outdir}/'3nv_Flagged_AllNonVariantsRepeatmaskQUAL_'${infile} \
-ef \
-o ${outdir}/'4nv_Filtered_AllNonVariantsRepeatmaskQUAL_'${infile}
## merge back together the variant and invariant files - hopefully we can do that.

java -jar -Xmx4G ${GATK} \
-T CombineVariants \
-R ${REFERENCE} \
-V ${outdir}/'4snp_VF_Filtered_GaTKHF_cluster_'${infile} \
-V ${outdir}/'4nv_Filtered_AllNonVariantsRepeatmaskQUAL_'${infile} \
-L ${mychrom} \
--assumeIdenticalSamples \
-o ${outdir}/'5_v4_mergeGaTKfiltered_varnonvar_'${infile}

## clean up
#rm ${outdir}/1_TrimAltRemoveNoCall_${infile}
#rm ${outdir}/'2snp_Filter_TrimAltRemoveNoCall_'${infile}
#rm ${outdir}/'3snp_VF_Flagged_GaTKHFGQ20_cluster_'${infile}
#rm ${outdir}/'4snp_VF_Filtered_GaTKHFGQ20_cluster_'${infile}
#rm ${outdir}/2nv_AllNonVariants_${infile}
#rm ${outdir}/'3nv_Flagged_AllNonVariantsRepeatmaskQUAL_'${infile}
#rm ${outdir}/'4nv_Filtered_AllNonVariantsRepeatmaskQUAL_'${infile}
