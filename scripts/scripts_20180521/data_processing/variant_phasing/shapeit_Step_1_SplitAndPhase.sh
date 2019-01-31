########## Shape it:

# install  v2.r904
#shapeit download (module on hoffman is too old)
#wget https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r904.glibcv2.12.linux.tar.gz
#tar -zxvf shapeit.v2.r904.glibcv2.12.linux.tar.gz

module load perl

shapeit=/u/home/a/ab08028/klohmueldata/annabel_data/bin/shapeit.v2.904.2.6.32-696.18.7.el6.x86_64/bin/shapeit
chromo=/u/home/a/ab08028/klohmueldata/annabel_data/bin/chromopainter-0.0.4/chromopainter

# added to path

wd=/u/flashscratch/a/ab08028/captures/
genotypedate=20181119
vcfdir=$wd/vcf_filtering/${genotypedate}_filtered
vcf=snp_7_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf.gz
phasedir=$wd/vcf_phasing/${genotypedate}_phased
mkdir -p $phasedir

# okay: missing data issue. either have to impute or remove missing data
# going to start by removing missing data and see what happens; downside: fewer snps , maybe conservative biased?
# upside: no imputation

######### NEED TO FIGURE OUT HOW TO SPLIT BY SCAFF: for now just testing all the way through with one scaff #########

# exclude relatives:
exclude=/u/flashscratch/a/ab08028/captures/samples/relativesToExcludeWhenPhasing.txt
# could also add high missingness exclude list 



############## check input #########
$shapeit check --input-vcf $vcfdir/$vcf --output-log $phasedir/${vcf%.vcf.gz}.stats --exclude-ind $exclude


# may need to deal with missingness inds
# can use   --exclude-ind
# note that two admixed inds have a lot of missing data (look at ind.mm file)

############ phase with shapeit  #############
$shapeit --thread 8 --input-vcf $vcfdir/$vcf -O $phasedir/${vcf%.vcf.gz} --window 0.5 --effective-size 4000 --exclude-ind $exclude
# 14 inds with high missingness >10%, leaving in with force for now
# 29000 snps included
# 4932 SNPs with high rates of missing data (>5%)
# only 29000 snps??
# 0.5 Mb window size
# pop effective size 4000 (should I combine across populations?)
# --rho constant recomb rate

# what does shapeit do to missing data?
############ convert to chromopainter format ############
# for now with backward compatibiltiy -f -v1 to work with chromopainter v1
perl impute2chromopainter.pl -f -v1 $phasedir/${vcf%.vcf.gz}.out.haps $phasedir/${vcf%.vcf.gz}.chromo.v1Format # seems to work (?)
# test with one scaffold  
# $shapeit --input-vcf test.vcf -O $phasedir/test.out --window 0.5 --effective-size 4000 --force 
# perl impute2chromopainter.pl test.out.haps chromoFormat.out # seems to work (?)

############### try chromo painter ############
###### stopped here -- need chromo2/globetrotter from garrett ######

################## optional convert output to vcf ########
# $shapeit -convert \
#--input-haps ${haplotypeOutFile} \
#--output-vcf ${vcfOutFile}
       
#issues: 
#threading
# split by chr

#######tests :

