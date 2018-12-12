#! /bin/bash
#$ -cwd
#$ -l h_rt=24:00:00,h_data=8G
#$ -N generate_cds_sfses
#$ -o /u/flashscratch/a/ab08028/captures/reports/SFS
#$ -e /u/flashscratch/a/ab08028/captures/reports/SFS
#$ -m abe
#$ -M ab08028
########### NOTE THAT VCF FILES MUST BE ANNOTATED USING VEP, pre-filtered for syn/miss using filter_vep and BGZIPPED #############
# this occurs in step flitering_Step_1_f-iii

# wrapper for making syn and missense SFSs with my custom script generate.StepByStep.SFS.FilterNoCallSimultaneously.py
# which for each site also removes the 
source /u/local/Modules/default/init/modules.sh
module load python/2.7
# sfs generating scripts:
gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptdir=$gitdir/scripts/scripts_20180521/analyses/generate_sfs/
script=generate.1D.SFS.FilterNoCallSimultaneously.py # this script will exclude any line that has no-call genotypes
# is different from regular script that will error out if there is a nocall genotype present because it indicates


genotypedate=20180806
vcfdir=/u/flashscratch/a/ab08028/captures/vcf_filtering/${genotypedate}_filtered/populationVCFs/cdsVCFs/
suffix='all_9_rmAllHet_rmRelativesAdmixed_passingAllFilters_allCalled'
prefix="all_9" # for output
populations="CA AK AL COM KUR"
outdir=/u/flashscratch/a/ab08028/captures/analyses/SFS/20180806/cdsSFS
mkdir -p $outdir

for type in syn missense
do
for pop in $populations
do

echo $pop
vcf=${pop}_VEP_${type}_${suffix}.vcf.gz

python $scriptdir/$script --vcf $vcfdir/$vcf --pop $pop --outdir $outdir --outPREFIX ${prefix}_${type}

done
done

