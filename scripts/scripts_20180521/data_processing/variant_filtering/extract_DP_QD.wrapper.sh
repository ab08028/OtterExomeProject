gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptdir=$gitdir/scripts/scripts_20180521/data_processing/variant_filtering

script=extract_DP.py 

genotypedate=20180806
wd=$SCRATCH/captures/vcf_filtering
vcfdir=$wd/${genotypedate}_filtered # date you called genotypes

vcf=all_1_TrimAlt_raw_variants.vcf.gz
scaff=GL896899.1
output=$vcfdir/${scaff}.${vcf%.vcf.gz}.DP.dist.txt
python $scriptdir/$script --VCF $vcfdir/$vcf --scaffold $scaff --outfile $output
gzip $output