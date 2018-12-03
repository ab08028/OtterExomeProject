gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptdir=$gitdir/scripts/scripts_20180521/data_processing/variant_filtering


genotypedate=20181119 # new 
wd=$SCRATCH/captures/vcf_filtering
vcfdir=$wd/${genotypedate}_filtered # date you called genotypes
rawvcfdir=/u/flashscratch/a/ab08028/captures/vcfs/vcf_20181119/

scaff=GL896899.1
mkdir -p $vcfdir/filteringStats
# get DP dist before any filtering (from raw file):
vcf1=raw_variants.vcf.gz

python $scriptdir/extract_DP.py  --VCF $rawvcfdir/$vcf1 --scaffold $scaff --outfile $vcfdir/filteringStats/${scaff}.${vcf1%.vcf.gz}.DP.dist.txt
gzip $vcfdir/filteringStats/${scaff}.${vcf1%.vcf.gz}.DP.dist.txt
# get QD dist from raw vcf:
python $scriptdir/extract_QD.py  --VCF $rawvcfdir/$vcf1 --scaffold $scaff --outfile $vcfdir/filteringStats/${scaff}.${vcf1%.vcf.gz}.QD.dist.txt
gzip $vcfdir/filteringStats/${scaff}.${vcf1%.vcf.gz}.QD.dist.txt


# get QD dist after min dp 500 filtering: 
vcf2=all_1_TrimAlt_raw_variants.vcf.gz

python $scriptdir/extract_QD.py  --VCF $vcfdir/$vcf2 --scaffold $scaff --outfile $vcfdir/filteringStats/${scaff}.${vcf2%.vcf.gz}.QD.dist.txt
gzip $vcfdir/${scaff}.${vcf2%.vcf.gz}.QD.dist.txt

# get QD for all snps across all scaffolds:
vcf3=snp_2_Filter_TrimAlt_raw_variants.vcf.gz 

python $scriptdir/extract_QD.allscaffs.py --VCF $vcfdir/$vcf3 --outfile $vcfdir/filteringStats/allscaffs.${vcf3%.vcf.gz}.QD.dist.txt


########### new script: pulls out QD, DP and QUAL for each site. ########
# get  dist from raw vcf:
python $scriptdir/extract_QD_DP_QUAL.py  --VCF $rawvcfdir/$vcf1 --scaffold $scaff --outfile $vcfdir/filteringStats/${scaff}.${vcf1%.vcf.gz}.joint.QD.DP.QUAL.dist.txt

# get dist from min dp 500 vcf:
python $scriptdir/extract_QD_DP_QUAL.py  --VCF $vcfdir/$vcf2 --scaffold $scaff --outfile $vcfdir/filteringStats/${scaff}.${vcf2%.vcf.gz}.joint.QD.DP.QUAL.dist.txt
