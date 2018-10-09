#! /bin/bash
#$ -cwd
#$ -l h_rt=24:00:00,h_data=8G,highp
#$ -N calculate_pi
#$ -o /u/flashscratch/a/ab08028/captures/reports/PI
#$ -e /u/flashscratch/a/ab08028/captures/reports/PI
#$ -m abe
#$ -M ab08028

### This is for the neutral SFS, but can modify choice of bed file to make coding SFS

source /u/local/Modules/default/init/modules.sh
module load python/2.7

rundate=20180806
### wrapper for Tanya's script
SCRATCH=/u/flashscratch/a/ab08028

# location of per-population input VCFs (with NO no-call genotypes)
vcfdir=$SCRATCH/captures/vcf_filtering/${rundate}_filtered

# temporary old neutral bed (probably can get more regions under new filtering scheme; this is a good conservative set)

neutralBed=${vcfdir}/bedCoords/all_8_rmRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_0.9_rmBadIndividuals_passingFilters.min10kb.fromExon.noCpGIsland.noRepeat.noFish.0based.sorted.merged.useThis.bed
mkdir -p ${vcfdir}/bedCoords/neutralCallableSites_perPop
# output SFS location
pidir=$SCRATCH/captures/analyses/pi/${rundate}
mkdir -p $pidir/neutralPi

# location of tanya's scripts
# this is latest 20180822 where it runs faster (but need to give enough memory)
# and doesn't require pre-filtered SFS 
tanyaDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/scripts_from_others/tanya_scripts/

# generate folded SFS:
populations="CA AK AL COM KUR"

for pop in $populations
do
echo $pop

# if vcf file is not prefiltered to just contain neutral (or other) regions
inVCF=${pop}_all_9_rmAllHet_rmRelativesAdmixed_passingAllFilters_allCalled.vcf.gz
# want to get callable sites in neutral regions per population; want a to be the vcf file because want to write out 
# this will find every site that is in the bed neutral regions that is present in the VCF
# it will output them as single sites in the bed:
# like this: GL896898.1	752539	752540
# and you can then sort and merge that bed to get new neutral regions for the specific pop
bedtools intersect -a $neutralBed -b $vcfdir/populationVCFs/$inVCF | bedtools sort -i stdin | bedtools merge -i stdin > ${vcfdir}/bedCoords/neutralCallableSites_perPop/${inVCF%.vcf.gz}.neutral.callableSites.0based.bed
popNeutBed=${vcfdir}/bedCoords/neutralCallableSites_perPop/${inVCF%.vcf.gz}.neutral.callableSites.0based.bed
python $tanyaDir/popgen_tools/popgen_tools.py \
--vcf_file $vcfdir/populationVCFs/$inVCF \
--pi_out $pidir/neutralPi/${inVCF%.vcf.gz}.pi.perPopCallableSiteNeutralRegions.out \
--total_SNPs $pidir/neutralPi/${inVCF%.vcf.gz}.totalSNPs.perPopCallableSiteNeutralRegions.out \
--no_sfs \
--target_bed $popNeutBed # instead of overall all-pop neutral regions, some of which may not be called int his population, use the actual ALL-CALLED neutral regions. will make pi calculations easier
done

# Note that now the output will be 
# start_region stop_region pi_in_region (not pi per site)
# note that if you use a general bed of neutral regions then this output won't represent callable sites
# but in this new way I'm doing it where I get the callable sites specifically beforehand, then it actually will
# give you the callable neutral sites that have passed all filters! 