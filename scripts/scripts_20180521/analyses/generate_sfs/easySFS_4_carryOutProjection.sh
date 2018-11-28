#! /bin/bash
#$ -cwd
#$ -l h_rt=30:00:00,h_data=32G,highp
#$ -N easySFSProjection
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
todaysdate=`date +%Y%m%d`
genotypeDate=20180806
vcfdir=/u/flashscratch/a/ab08028/captures/vcf_filtering/${genotypeDate}_filtered/neutralVCF_allPops
popFile=/u/flashscratch/a/ab08028/captures/samples/samplesPop.Headers.forEasySFS.2.txt
# this has admixed in it , but they aren't in pop file
#easySFS=/u/home/a/ab08028/klohmueldata/annabel_data/bin/easySFS/easySFS.abContinueMod.py 



easySFS=/u/home/a/ab08028/klohmueldata/annabel_data/bin/easySFS/easySFS.abModified.2.noInteract.Exclude01Sites.20181121.py  # this is my modification
# this version of script excludes sites that are 0-1 across all populations (maybe) -- not sure if it does yet. 

## choose your projections:
projections="20,20,20,20,20"

outdir=/u/flashscratch/a/ab08028/captures/analyses/SFS/$genotypeDate/easySFS
mkdir -p $outdir
# had to modify easySFS so that it wouldn't prompt a "yes/no" response about samples that are missing from VCF file

# make sure vcf isn't zipped

#vcf=neutral.snp_8a_rmRelatives_rmAdmixedOutliers_passingBespoke_maxNoCallFrac_0.9_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf
vcf='neutral.all_8a_rmRelatives_rmAdmixedOutliers_passingBespoke_maxNoCallFrac_'${noCallFrac}'_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf'
# hoping I can get it to work with this neutral all-sites file so that I have monomorphic sites sampled as well. This would be ideal.
gunzip ${vcf}.gz
$easySFS -i $vcfdir/${vcf} -p $popFile -a -v --proj 20,20,20,20,20 -f -o $outdir/projection-${todaysdate}
# test run only had --proj 20 (maybe just projects 1 population? will have to see)
# -f forces overwrite of outdir
$bgzip ${vcf}
# then do for SYN and MIS (eventually)
