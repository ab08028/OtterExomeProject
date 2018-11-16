####### IGV wrapper:
source /u/local/Modules/default/init/modules.sh

module load python/2.7
module load samtools

todaysdate=2018113
gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptdir=$gitdir/scripts/scripts_20180521/data_processing/variant_filtering/sandbox
script=getRegionsForIGV.py
bamdir=/u/flashscratch/a/ab08028/captures/paleomix/
igvdir=/u/flashscratch/a/ab08028/captures/IGV/

######### files ########
rundate=20180806 # date genotypes were called (vcf_20180806 includes capture 02)
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/vcf_filtering
vcfdir=$wd/${rundate}_filtered/populationVCFs/neutralVCFs # date you called genotypes



vcf=$vcfdir/CA_neutral_all_9_rmAllHet_rmRelativesAdmixed_passingAllFilters_allCalled.vcf.gz
 # temp
# usage
# python getRegionsForIGV.py [infile full path] [outfile script path] [path to paleomix dir where bams are] [path to where you want new bams to go] 

python $scriptdir/$script $vcf $scriptdir/gatherRegionsForIGV.${todaysdate}.sh $bamdir $igvdir
# this generates the script

# then run the script:
sh $scriptdir/gatherRegionsForIGV.${todaysdate}.sh

# need to index the bams