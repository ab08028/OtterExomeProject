#! /bin/bash
#$ -cwd
#$ -l h_rt=24:00:00,h_data=8G
#$ -N dadi_inference
#$ -o /u/flashscratch/a/ab08028/captures/reports/cdsSites
#$ -e /u/flashscratch/a/ab08028/captures/reports/cdsSites
#$ -m abe
#$ -M ab08028

### count callable sites per individual:
# already counted mis/syn per individual
# just need callable cds:
source /u/local/Modules/default/init/modules.sh
module load python
######### directories and files #########
SCRATCH=/u/flashscratch/a/ab08028/
genotypedate=20181119
vcfdir=$SCRATCH/captures/vcf_filtering/${genotypedate}_filtered/neutral_and_cds_VCFs/cdsVCFs
gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject
scriptdir=$gitdir/scripts/scripts_20180521/analyses/compareMissense_Syn
script=countCallableSitesPerIndividual.py
allSitesCDSVCF=cds_all_9_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters_raw_variants.vcf.gz
# this vcf is CDS Only(!!); admixed and rels removed, and 75% het filter applied on top of all other filters. No missing data filter applied, however (may bias things? we shall see)
outdir=/u/flashscratch/a/ab08028/captures/analyses/compareMissense_Syn/callableCDSPerIndividual
mkdir -p $outdir
todaysdate=`date +%Y%m%d`

### run script: 
python $scriptdir/$script --vcf $vcfdir/$allSitesCDSVCF --outfile $outdir/callableCDSSitesPerIndividual.${todaysdate}.txt

# should output a table with called and missing gts per ind (called+missing should be the same for all)
