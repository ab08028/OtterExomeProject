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
tanyaDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/data_processing/tanya_scripts/

# neutral sites that have been called (min 10kb from genes)
neutralBed=${vcfDir}/bedCoords/all_7_passingBespoke.min10kb.fromExon.0based.sorted.merged.bed

# generate folded SFS:
populations="CA AK AL COM KUR"

for pop in $populations
do
inVCF=${pop}_all_7_passingAllFilters_allCalledraw_variants.vcf.gz

python $tanyaDir/popgen_tools/popgen_tools.py \
--vcf_file $vcfdir/populationVCFs/$inVCF \
--target_bed $neutralBed \
--total_SNPs $SFSdir/${inVCF%.vcf.gz}.sfs.out \
--no_pi

done