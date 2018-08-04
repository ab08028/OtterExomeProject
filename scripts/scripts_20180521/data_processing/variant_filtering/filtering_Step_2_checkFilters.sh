# check that alt and ref are not indels

rundate=20180724 # date genotypes were called

#### file locations
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/vcf_filtering
mkdir -p $wd
indir=$SCRATCH/captures/vcfs/vcf_${rundate}
infile=raw_variants.vcf.gz ### make sure this doesn't have a path as part of its name! just infile names
outdir=$indir/checkFilters
mkdir -p $outdir 
zcat all_5_passingFilters_80percCall_${infile} | grep -v "#" | awk 'length($4)>1 || length($5)>1' > indelLines.all.test.vcf
zcat snp_5_passingAllFilters_postMerge_raw_variants.${infile} | grep -v "#" | awk 'length($4)>1 || length($5)>1' > indelLines.snps.vcf

# are all hets: grep "0/1" gets your hets (****assumes that data is unphased **** code this more defensively )
zcat all_5_passingFilters_80percCall_raw_variants.vcf.gz | grep "0/1" | grep -v "1/1" | grep -v "0/0" | wc -l
