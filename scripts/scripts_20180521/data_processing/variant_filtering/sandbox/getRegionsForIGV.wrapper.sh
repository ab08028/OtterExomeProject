####### IGV wrapper:

module load python/2.7
todaysdate=2018113
gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptdir=$gitdir/scripts/scripts_20180521/data_processing/variant_filtering/
script=getRegionsForIGV.py
bamdir=/u/flashscratch/a/ab08028/captures/paleomix/
igvdir=/u/flashscratch/a/ab08028/captures/IGV/
# usage
# python getRegionsForIGV.py [infile full path] [outfile script path] [path to paleomix dir where bams are] [path to where you want new bams to go] 

python $script $vcf $scriptdir/gatherRegionsForIGV.${todaysdate}.sh $bamdir $igvdir