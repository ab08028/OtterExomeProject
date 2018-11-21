#! /bin/bash
#$ -cwd
#$ -l h_rt=30:00:00,h_data=32G,highp
#$ -N easySFSPreview
#$ -o /u/flashscratch/a/ab08028/captures/reports/SFS
#$ -e /u/flashscratch/a/ab08028/captures/reports/SFS
#$ -m abe
#$ -M ab08028

####### Easy SFS
# https://github.com/isaacovercast/easySFS
# install:
# git clone git@github.com:isaacovercast/easySFS.git
# cd easySFS
# chmod +x *.py
# easySFS.py
source /u/local/Modules/default/init/modules.sh
module load python/2.7
bgzip=/u/home/a/ab08028/klohmueldata/annabel_data/bin/tabix-0.2.6/bgzip

genotypeDate=20180806
vcfdir=/u/flashscratch/a/ab08028/captures/vcf_filtering/${genotypeDate}_filtered/
#vcf=all_8_rmRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_0.9_rmBadIndividuals_passingFilters_raw_variants.vcf.gz
popFile=/u/flashscratch/a/ab08028/captures/samples/samplesPop.Headers.forEasySFS.2.txt
#vcf=snp_8_rmRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf # had to unzip it
#vcfgz=all_8_rmRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_0.9_rmBadIndividuals_passingFilters_raw_variants.vcf.gz
#gunzip $vcfgz
#vcf=${vcfgz%.gz}
# this has admixed in it , but they aren't in pop file
easySFS=/u/home/a/ab08028/klohmueldata/annabel_data/bin/easySFS/easySFS.abContinueMod.py
outdir=/u/flashscratch/a/ab08028/captures/analyses/SFS/$genotypeDate/easySFS
mkdir -p $outdir
# had to modify easySFS so that it wouldn't prompt a "yes/no" response about samples that are missing from VCF file
# want it to just output that info and continue , not prompt yes/no.
# test:

# need to get it working for multi sfs 

#populations="CA AK AL COM KUR"
#for pop in $populations
#do
# this vcf file has all snps with maxnocallfrac 0.9 (liberal), and with admixed and outliers and relatives all removed
# (so it's the samples that will go into the SFS)

vcf=snp_8a_rmRelatives_rmAdmixedOutliers_passingBespoke_maxNoCallFrac_0.9_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf
#gunzip $vcfdir/populationVCFs/$vcf it must be unzipped
#( you are here )
$easySFS -i $vcfdir/${vcf} -p $popFile --preview -a -v > $outdir/${vcf%.vcf}.easySFS.projPreview.txt

# $bgzip $vcfdir/${vcf%.gz}
# done

# need to add -a otherwise it selects one snp per locus (was built for RAD data)
# eventually make popFile list consistent with final dataset individuals
# otherwise have to say "yes" on screen
    #Running preview mode. We will print out the results for # of segregating sites
    #for multiple values of projecting down for each population. The dadi
    #manual recommends maximizing the # of seg sites for projections, but also
   	#a balance must be struck between # of seg sites and sample size.

    #For each population you should choose the value of the projection that looks
   # best and then rerun easySFS with the `--proj` flag.
# test: needs to not be compressed?
###### Runs out of memory!!! Try on SNP file instead, since want to maximize SNPs anyway. but snp file already has a 0.2 filter applied
# which we don't want 
