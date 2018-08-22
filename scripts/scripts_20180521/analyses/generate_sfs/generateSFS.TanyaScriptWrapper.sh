#! /bin/bash
#$ -cwd
#$ -l h_rt=50:00:00,h_data=16G,highp
#$ -N generate_sfs_allPops
#$ -o /u/flashscratch/a/ab08028/captures/reports/SFS
#$ -e /u/flashscratch/a/ab08028/captures/reports/SFS
#$ -m abe
#$ -M ab08028


source /u/local/Modules/default/init/modules.sh
module load python/2.7

rundate=20180806
### wrapper for Tanya's script
SCRATCH=/u/flashscratch/a/ab08028

# location of per-population input VCFs (with NO no-call genotypes)
vcfdir=$SCRATCH/captures/vcf_filtering/${rundate}_filtered

# output SFS location
SFSdir=$SCRATCH/captures/analyses/SFS/${rundate}
mkdir -p $SFSdir

# location of tanya's scripts
tanyaDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/scripts_from_others/tanya_scripts/

# generate folded SFS:
populations="CA AK AL COM KUR"

for pop in $populations
do
echo $pop
inVCF=${pop}_neutral_7_passingAllFilters_allCalledraw_variants.vcf.gz # prefiltered to only be neutral regions!

# build SFS
python $tanyaDir/popgen_tools/popgen_tools.py \
--vcf_file $vcfdir/populationVCFs/$inVCF \
--total_SNPs $SFSdir/${inVCF%.vcf.gz}.sfs.out \
--no_pi \
--sfs_no_target_bed 
# added --sfs_no_target_bed so that it doesn't cross ref every region with bed file
# instead, pre filter the vcf to only be neutral regions

done
