#! /bin/bash
#$ -cwd
#$ -l h_rt=24:00:00,h_data=8G
#$ -N generate_sfs_allPops_afterFiltering
#$ -o /u/flashscratch/a/ab08028/captures/reports/SFS
#$ -e /u/flashscratch/a/ab08028/captures/reports/SFS
#$ -m abe
#$ -M ab08028

# wrapper for making neutral SFSs with my custom script generate1DSFS.py
source /u/local/Modules/default/init/modules.sh
module load python/2.7
gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptdir=$gitdir/scripts/scripts_20180521/analyses/generate_sfs/
script=generate1DSFS.py

genotypedate=20180806
vcfdir=/u/flashscratch/a/ab08028/captures/vcf_filtering/${genotypedate}_filtered/populationVCFs/neutralVCFs/
suffix='neutral_all_9_rmAllHet_rmRelativesAdmixed_passingAllFilters_allCalled.vcf.gz'
populations="CA AK AL COM KUR"
outdir=/u/flashscratch/a/ab08028/captures/analyses/SFS/20180806/neutralSFS

for pop in $populations
do
echo $pop
python $scriptdir/$script --vcf $vcfdir/${pop}_${suffix} --pop $pop --outdir $outdir --outPREFIX "all_9"
done

# skipping admixed SFSs for now
