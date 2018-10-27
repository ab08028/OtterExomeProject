###### vcftools pi calculations (can be run in the shell or as a job)
source /u/local/Modules/default/init/modules.sh
module load vcftools
######################## neutral pi #######################
# calculate it based on the neutral VCF files
genotypeDate=20180806
wd=/u/flashscratch/a/ab08028/captures/
vcfdir=$wd/vcf_filtering/${genotypeDate}_filtered/populationVCFs/neutralVCFs
outdir=$wd/analyses/pi/${genotypeDate}/neutralPi

suffix=neutral_all_9_rmAllHet_rmRelativesAdmixed_passingAllFilters_allCalled
for pop in CA AK AL COM KUR
do
vcf=${pop}_${suffix}.vcf.gz
zcat $vcfdir/$vcf | vcftools --vcf - --site-pi  --out ${outdir}/${pop}_${suffix}
# note "-" after --vcf
gzip ${outdir}/${pop}_${suffix}.sites.pi
done




 