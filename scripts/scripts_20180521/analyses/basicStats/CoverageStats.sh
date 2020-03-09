#$ -N gatherDepthStats
#$ -o /u/flashscratch/a/ab08028/captures/reports/depth
#$ -e /u/flashscratch/a/ab08028/captures/reports/depth
#$ -m abe
#$ -M ab08028
#$ -l h_rt=24:00:00,h_data=20G

######## get some depth stats
# want to do this on all sites (not snps)
#filtering
module load vcftools
vcfdir=/u/flashscratch/a/ab08028/captures/vcf_filtering/20181119_filtered
outdir=$vcfdir/depthAndMissingnessStats_20200309 
############## start with all_7 ##################
vcf=all_7_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters_raw_variants.vcf.gz
prefix=all_7
vcftools --gzvcf $vcfdir/$vcf --out $outdir/$prefix.vcftoolsStats --depth --missing-indv
mv $outdir/$prefix.vcftoolsStats $outdir/$prefix.vcftoolsStats.imiss $outdir/$prefix.vcftoolsStats.MissingnessPerIndividual.imiss.txt
mv $outdir/$prefix.vcftoolsStats $outdir/$prefix.vcftoolsStats.idepth $outdir/$prefix.vcftoolsStats.MeanDepthPerIndividual.idepth.txt
zcat $vcfdir/$vcf | grep -v "#" -c > $outdir/$prefix.TOTALSITESINVCF.txt

##### skip: all_8:  didn't change filters, just removed relatives and admixed (so shouldn't be needed) #######
#vcf=all_8_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters_raw_variants.vcf.gz 

###### all_9: added het 0.75 filter  ; no missingness filter #######
vcf=all_9_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters_raw_variants.vcf.gz
prefix=all_9
vcftools --gzvcf $vcfdir/$vcf --out $outdir/$prefix.ccftoolsStats --depth --missing-indv
mv $outdir/$prefix.vcftoolsStats $outdir/$prefix.vcftoolsStats.imiss $outdir/$prefix.vcftoolsStats.MissingnessPerIndividual.imiss.txt
mv $outdir/$prefix.vcftoolsStats $outdir/$prefix.vcftoolsStats.idepth $outdir/$prefix.vcftoolsStats.MeanDepthPerIndividual.idepth.txt
# count up total sites
zcat $vcfdir/$vcf | grep -v "#" -c > $outdir/$prefix.TOTALSITESINVCF.txt