#! /bin/bash
#$ -cwd
#$ -l h_rt=24:00:00,h_data=16G
#$ -N vcf1c_filtering
#$ -o /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -e /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -m abe
#$ -M ab08028

# modules
source /u/local/Modules/default/init/modules.sh
module load java
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

# location of vcf checking and filtering script
# incompatible scaffolds: repeatMaskCoords=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/repeatMaskingCoordinates/masking_coordinates.bed
# repeat masking is optional: my target captrue was designed away from repeat regions, so not a huge deal 
# if you don't include this; but trying to be extra thorough.
# these coords are downloaded from NCBI FTP site for mustela putorius furo ; if using elut genome, use Annotation repeat coords; 
# ftp://ftp.ncbi.nlm.nih.gov/genomes/Mustela_putorius_furo/ <-- location of mustela files
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Mustela_putorius_furo/masking_coordinates.gz

outdir=$wd/${rundate}_filtered # date you called genotypes
mkdir -p $outdir
##################################################################
############################ GET STATS ###########################
##################################################################
mkdir -p ${outdir}/filteringStats/
echo "starting getting statistics"
# starting sites:
echo "# Filtering Statistics for SNPs and invariant sites" > ${outdir}/filteringStats/filteringStats.${rundate}.txt
stat0=`zcat ${indir}/${infile} | grep -v -c "#"` 
echo "stat0 starting_sites" $stat0 >> ${outdir}/filteringStats/filteringStats.${rundate}.txt
# sites after trim Altnernates and removing sites with >20% no call
stat1=`zcat ${outdir}/all_1_TrimAlt_${infile} | grep -v -c "#"`
echo "stat1 sites_afterTrimAlt" $stat1 >> ${outdir}/filteringStats/filteringStats.${rundate}.txt

########## biallelic snps: ##########
stat2=`zcat ${outdir}/'snp_2_Filter_TrimAlt_'${infile} | grep -v -c "#"`
echo "stat2 biSNPS_initial" $stat2 >> ${outdir}/filteringStats/filteringStats.${rundate}.txt

# snps failing filters:
# note these may overlap (same snp can fail multiple things)

stat31=`zcat ${outdir}/'snp_3_Flagged_GQ_DP_GaTKHF_cluster_'${infile} | grep -v "#" | grep -c FAIL_RepMask`
echo "stat3.1 SNPS_failingREPMASK" $stat31 >> ${outdir}/filteringStats/filteringStats.${rundate}.txt
stat32=`zcat ${outdir}/'snp_3_Flagged_GQ_DP_GaTKHF_cluster_'${infile} | grep -v "#" | grep -c FAIL_GATKHF`
echo "stat3.2 SNPS_failingGATKHF" $stat32 >> ${outdir}/filteringStats/filteringStats.${rundate}.txt
stat36=`zcat ${outdir}/'snp_3_Flagged_GQ_DP_GaTKHF_cluster_'${infile} | grep -v "#" | grep -c SnpCluster`
echo "stat3.6 SNPS_failingSNPCluster" $stat36 >> ${outdir}/filteringStats/filteringStats.${rundate}.txt
# genotypes that fail:
# note that grep -o outputs ALL occurences, not all lines (so you can count multiple on one line)
stat33=`zcat ${outdir}/'snp_3_Flagged_GQ_DP_GaTKHF_cluster_'${infile} | grep -v "#" | grep -o FAIL_GQ | wc -l`
echo "stat3.3 var_GTs_failingGQ" $stat33 >> ${outdir}/filteringStats/filteringStats.${rundate}.txt
stat34=`zcat ${outdir}/'snp_3_Flagged_GQ_DP_GaTKHF_cluster_'${infile} | grep -v "#" | grep -o FAIL_DP_HIGH | wc -l`
echo "stat3.4 var_GTs_failingDP_HIGH" $stat34 >> ${outdir}/filteringStats/filteringStats.${rundate}.txt
stat35=`zcat ${outdir}/'snp_3_Flagged_GQ_DP_GaTKHF_cluster_'${infile} | grep -v "#" | grep -o FAIL_DP_LOW| wc -l`
echo "stat3.5 var_GTs_failingDP_LOW" $stat35 >> ${outdir}/filteringStats/filteringStats.${rundate}.txt

# total passing snps
stat4=`zcat ${outdir}/'snp_4_Filtered_GQ_DP_GaTKHF_cluster_'${infile} | grep -v -c "#"`
echo "stat4 totalPassingSNPsSites_preMerge" $stat4 >> ${outdir}/filteringStats/filteringStats.${rundate}.txt

########## invariant sites: ###########
stat5=`zcat ${outdir}/'nv_2_AllNonVariants_'${infile} | grep -v -c "#"` 
echo "stat5 invarSites_initial" $stat5 >> ${outdir}/filteringStats/filteringStats.${rundate}.txt

# filters of invariant sites:
stat61=`zcat ${outdir}/'nv_3_Flagged_DP_RGQ_QUAL_'${infile} | grep -v "#" | grep -c FAIL_RepMask`
echo "stat6.1 invarSites_failingREPMASK" $stat61 >> ${outdir}/filteringStats/filteringStats.${rundate}.txt
stat62=`zcat ${outdir}/'nv_3_Flagged_DP_RGQ_QUAL_'${infile} | grep -v "#" | grep -c FAIL_QUAL30`
echo "stat6.2 invarSites_failingQUAL30" $stat62 >> ${outdir}/filteringStats/filteringStats.${rundate}.txt
# invariant genotype filters:
stat63=`zcat ${outdir}/'nv_3_Flagged_DP_RGQ_QUAL_'${infile} | grep -v "#" | grep -o FAIL_DP_LOW | wc -l`
echo "stat6.3 invar_GTs_failingDP_LOW" $stat63 >> ${outdir}/filteringStats/filteringStats.${rundate}.txt
stat64=`zcat ${outdir}/'nv_3_Flagged_DP_RGQ_QUAL_'${infile} | grep -v "#" | grep -o FAIL_DP_HIGH | wc -l`
echo "stat6.4 invar_GTs_failingDP_HIGH" $stat64 >> ${outdir}/filteringStats/filteringStats.${rundate}.txt
stat65=`zcat ${outdir}/'nv_3_Flagged_DP_RGQ_QUAL_'${infile} | grep -v "#" | grep -o FAIL_RGQ | wc -l`
echo "stat6.5 invar_GTs_failingRGQ" $stat65 >> ${outdir}/filteringStats/filteringStats.${rundate}.txt

# total passing invariant sites 
stat7=`zcat ${outdir}/'nv_4_Filtered_DP_RGQ_QUAL_'${infile} | grep -v -c "#"`
echo "stat7 total_invarSites_passing_preMerge" $stat7 >> ${outdir}/filteringStats/filteringStats.${rundate}.txt

################## merged sites ##############
# total passing sites
stat81=`zcat ${outdir}/'all_5_passingFilters_'${infile} | grep -v -c "#"`
echo "stat8.1 totalPassingSites_all_preBespoke" $stat81 >> ${outdir}/filteringStats/filteringStats.${rundate}.txt

stat82=`zcat ${outdir}/'all_7_passingBespoke_passingFilters_'${infile} | grep -v -c "#"`
echo "stat8.2 totalPassingSites_all_postBespoke" $stat82 >> ${outdir}/filteringStats/filteringStats.${rundate}.txt

stat83=`grep -v -c "#" ${outdir}/'fail_all_7_FAILINGBespoke_passingFilters_'${infile%.vcf.gz}.txt`
echo "stat8.3 sites_failing_bespoke" $stat83 >> ${outdir}/filteringStats/filteringStats.${rundate}.txt

# note that some snps may have become invariant after filtering.
stat9=`zcat ${outdir}/'snp_7_80perc_passingBespoke_passingAllFilters_postMerge_'${infile} | grep -v -c "#"` 
echo "stat9 totalPassingSNPsSites_80Perc_postMerge_postBespoke" $stat9 >> ${outdir}/filteringStats/filteringStats.${rundate}.txt

stat10=`zcat ${outdir}/'nv_7_80perc_passingBespoke_passingAllFilters_postMerge_'${infile} | grep -v -c "#"` 
echo "stat10 total_invarSites_passing_80Perc_postMerge_postBespoke" $stat10 >> ${outdir}/filteringStats/filteringStats.${rundate}.txt
echo "done getting statistics"