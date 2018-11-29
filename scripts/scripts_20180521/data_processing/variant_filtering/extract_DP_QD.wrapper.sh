gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptdir=$gitdir/scripts/scripts_20180521/data_processing/variant_filtering


genotypedate=20180806
wd=$SCRATCH/captures/vcf_filtering
vcfdir=$wd/${genotypedate}_filtered # date you called genotypes

vcf=all_1_TrimAlt_raw_variants.vcf.gz
scaff=GL896899.1

# get DP dist:

python $scriptdir/extract_DP.py  --VCF $vcfdir/$vcf --scaffold $scaff --outfile $vcfdir/${scaff}.${vcf%.vcf.gz}.DP.dist.txt
gzip $vcfdir/${scaff}.${vcf%.vcf.gz}.DP.dist.txt


# get QD dist: 
python $scriptdir/extract_QD.py  --VCF $vcfdir/$vcf --scaffold $scaff --outfile $vcfdir/${scaff}.${vcf%.vcf.gz}.QD.dist.txt
gzip $vcfdir/${scaff}.${vcf%.vcf.gz}.QD.dist.txt

# get QD dist POST min depth 500 filtering:
vcf2=snp_2_Filter_TrimAlt_raw_variants.vcf.gz
python $scriptdir/extract_QD.py  --VCF $vcfdir/$vcf2 --scaffold $scaff --outfile $vcfdir/${scaff}.${vcf%.vcf.gz}.QD.dist.txt
