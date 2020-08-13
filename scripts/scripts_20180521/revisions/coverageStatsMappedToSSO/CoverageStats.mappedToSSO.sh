#$ -N gatherDepthStats
#$ -o /u/flashscratch/a/ab08028/captures/reports/depth
#$ -e /u/flashscratch/a/ab08028/captures/reports/depth
#$ -m abe
#$ -M ab08028
#$ -l h_rt=24:00:00,h_data=20G

######## get some depth stats
# want to do this on all sites (not snps)
#filtering
source /u/local/Modules/default/init/modules.sh
module load vcftools
vcfdir=/u/flashscratch/a/ab08028/captures/vcf_filtering/20200719_SSO_filtered
outdir=$vcfdir/depthAndMissingnessStats
mkdir -p $outdir

############## start with all_7 ##################
vcf=all_7_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters_raw_variants.vcf.gz
prefix=all_7
vcftools --gzvcf $vcfdir/$vcf --out $outdir/$prefix.vcftoolsStats.MeanDepthPerIndividual --depth # get mean depth per ind 
zcat $vcfdir/$vcf | grep -v "#" -c > $outdir/$prefix.TOTALSITESINVCF.txt # get total sites 

