#! /bin/bash
#$ -cwd
#$ -l h_rt=24:00:00,h_data=16G
#$ -N het_per_individual
#$ -o /u/flashscratch/a/ab08028/captures/reports/het
#$ -e /u/flashscratch/a/ab08028/captures/reports/het
#$ -m abe
#$ -M ab08028
###### Het per individual wrapper
source /u/local/Modules/default/init/modules.sh
module load python/2.7
REFSHORTCODE=SSO # mapped to southern sea otter 
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/vcf_filtering
scriptdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_filtering
script=filtering_getHetPerInd.py
rundate=20200719_${REFSHORTCODE}
# for all sites:
#vcfdir=$wd/${rundate}_filtered
# for neutral:
vcfdir=$wd/${rundate}_filtered/neutral/neutral_and_cds_VCFs/neutralVCFs
outdir=$vcfdir/hetPerIndividual/$REFSHORTCODE
mkdir -p $outdir

# there are a two vcfs I want to look at and compare
#vcfs='all_7_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters_raw_variants.vcf.gz all_9_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters_raw_variants.vcf.gz'
#vcfs='all_9_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters_raw_variants.vcf.gz'
vcfs='neutral.all_9_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters_raw_variants.vcf.gz'
for vcf in $vcfs
do
python $scriptdir/$script $vcfdir/$vcf $outdir/${vcf}.perIndHet.${REFSHORTCODE}.txt
done
