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

# incompatible scaffolds: repeatMaskCoords=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/repeatMaskingCoordinates/masking_coordinates.bed
# repeat masking is optional: my target captrue was designed away from repeat regions, so not a huge deal 
# if you don't include this; but trying to be extra thorough.
# these coords are downloaded from NCBI FTP site for mustela putorius furo ; if using elut genome, use Annotation repeat coords; 
# ftp://ftp.ncbi.nlm.nih.gov/genomes/Mustela_putorius_furo/ <-- location of mustela files
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Mustela_putorius_furo/masking_coordinates.gz

outdir=$wd/${rundate}_filtered # date you called genotypes

########## biallelic snps: ##########
stat2=`zcat ${outdir}/'snp_2_Filter_TrimAlt80Perc_'${infile} | grep -v -c "#"`
echo "stat2 biSNPS_initial" $stat2 >> ${outdir}/filteringStats.${rundate}.txt

# snps failing filters:
# note these may overlap (same snp can fail multiple things)

#stat31=`zcat ${outdir}/'snp_3_Flagged_GQ_DP_GaTKHF_cluster_'${infile} | grep -v "#" | grep -c FAIL_RepMask`
#echo "stat3.1 SNPS_failingREPMASK" $stat31 >> ${outdir}/filteringStats.${rundate}.txt
#stat32=`zcat ${outdir}/'snp_3_Flagged_GQ_DP_GaTKHF_cluster_'${infile} | grep -v "#" | grep -c FAIL_GATKHF`
#echo "stat3.2 SNPS_failingGATKHF" $stat32 >> ${outdir}/filteringStats.${rundate}.txt
#stat36=`zcat ${outdir}/'snp_3_Flagged_GQ_DP_GaTKHF_cluster_'${infile} | grep -v "#" | grep -c SnpCluster`
#echo "stat3.6 SNPS_failingSNPCluster" $stat36 >> ${outdir}/filteringStats.${rundate}.txt
# genotypes that fail:
# note that grep -o outputs ALL occurences, not all lines (so you can count multiple on one line)
#stat33=`zcat ${outdir}/'snp_3_Flagged_GQ_DP_GaTKHF_cluster_'${infile} | grep -v "#" | grep -o FAIL_GQ | wc -l`
#echo "stat3.3 var_GTs_failingGQ" $stat33 >> ${outdir}/filteringStats.${rundate}.txt
#stat34=`zcat ${outdir}/'snp_3_Flagged_GQ_DP_GaTKHF_cluster_'${infile} | grep -v "#" | grep -o FAIL_DP_HIGH | wc -l`
#echo "stat3.4 var_GTs_failingDP_HIGH" $stat34 >> ${outdir}/filteringStats.${rundate}.txt
#stat35=`zcat ${outdir}/'snp_3_Flagged_GQ_DP_GaTKHF_cluster_'${infile} | grep -v "#" | grep -o FAIL_DP_LOW| wc -l`
#echo "stat3.5 var_GTs_failingDP_LOW" $stat35 >> ${outdir}/filteringStats.${rundate}.txt

# total passing snps
stat4=`zcat ${outdir}/'snp_4_Filtered_80percCall_GQ_DP_GaTKHF_cluster_'${infile} | grep -v -c "#"`
echo "stat4 totalPassingSNPsSites_80Perc_preMerge" $stat4 >> ${outdir}/filteringStats.${rundate}.txt

########## invariant sites: ###########
stat5=`zcat ${outdir}/'nv_2_AllNonVariants_'${infile} | grep -v -c "#"` 
echo "stat5 invarSites_initial" $stat5 >> ${outdir}/filteringStats.${rundate}.txt

# filters of invariant sites:
#stat61=`zcat ${outdir}/'nv_3_Flagged_DP_RGQ_QUAL_'${infile} | grep -v "#" | grep -c FAIL_RepMask`
#echo "stat6.1 invarSites_failingREPMASK" $stat61 >> ${outdir}/filteringStats.${rundate}.txt
#stat62=`zcat ${outdir}/'nv_3_Flagged_DP_RGQ_QUAL_'${infile} | grep -v "#" | grep -c FAIL_QUAL30`
#echo "stat6.2 invarSites_failingQUAL30" $stat62 >> ${outdir}/filteringStats.${rundate}.txt
# invariant genotype filters:
#stat63=`zcat ${outdir}/'nv_3_Flagged_DP_RGQ_QUAL_'${infile} | grep -v "#" | grep -o FAIL_DP_LOW | wc -l`
#echo "stat6.3 invar_GTs_failingDP_LOW" $stat63 >> ${outdir}/filteringStats.${rundate}.txt
#stat64=`zcat ${outdir}/'nv_3_Flagged_DP_RGQ_QUAL_'${infile} | grep -v "#" | grep -o FAIL_DP_HIGH | wc -l`
#echo "stat6.4 invar_GTs_failingDP_HIGH" $stat64 >> ${outdir}/filteringStats.${rundate}.txt
#stat65=`zcat ${outdir}/'nv_3_Flagged_DP_RGQ_QUAL_'${infile} | grep -v "#" | grep -o FAIL_RGQ | wc -l`
#echo "stat6.5 invar_GTs_failingRGQ" $stat65 >> ${outdir}/filteringStats.${rundate}.txt

# total passing invariant sites 
stat7=`zcat ${outdir}/'nv_4_Filtered_80percCall_DP_RGQ_QUAL_'${infile} | grep -v -c "#"`
echo "stat7 total_invarSites_passing_80Perc_preMerge" $stat7 >> ${outdir}/filteringStats.${rundate}.txt

################## merged sites ##############
# total passing sites
stat8=`zcat ${outdir}/'all_5_passingFilters_80percCall_'${infile} | grep -v -c "#"`
echo "stat8 totalPassingSites_all" $stat8 >> ${outdir}/filteringStats.${rundate}.txt
# note that some snps may have become invariant after filtering.
stat9=`zcat ${outdir}/'snp_5_passingAllFilters_postMerge_'${infile} | grep -v -c "#"` 
echo "stat9 totalPassingSNPsSites_80Perc_postMerge" $stat9 >> ${outdir}/filteringStats.${rundate}.txt

stat10=`zcat ${outdir}/'nv_5_passingAllFilters_postMerge_'${infile} | grep -v -c "#"` 
echo "stat10 total_invarSites_passing_80Perc_postMerge" $stat10 >> ${outdir}/filteringStats.${rundate}.txt
echo "done getting statistics"

## optional : clean up
#rm ${outdir}/1_TrimAltRemoveNoCall_${infile}
#rm ${outdir}/'2snp_Filter_TrimAltRemoveNoCall_'${infile}
#rm ${outdir}/'3snp_VF_Flagged_GaTKHFGQ20_cluster_'${infile}
#rm ${outdir}/'4snp_VF_Filtered_GaTKHFGQ20_cluster_'${infile}
#rm ${outdir}/2nv_AllNonVariants_${infile}
#rm ${outdir}/'3nv_Flagged_AllNonVariantsRepeatmaskQUAL_'${infile}
#rm ${outdir}/'4nv_Filtered_AllNonVariantsRepeatmaskQUAL_'${infile}


