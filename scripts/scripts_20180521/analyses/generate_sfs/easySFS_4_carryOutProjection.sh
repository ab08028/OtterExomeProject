#! /bin/bash
#$ -cwd
#$ -l h_rt=30:00:00,h_data=64G,highp
#$ -N easySFSProjection2
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
todaysdate=20181127-64gb
genotypeDate=20180806
noCallFrac=0.9
vcfdir=/u/flashscratch/a/ab08028/captures/vcf_filtering/${genotypeDate}_filtered/neutralVCF_allPops
popFile=/u/flashscratch/a/ab08028/captures/samples/samplesPop.Headers.forEasySFS.2.txt
# this has admixed in it , but they aren't in pop file
#easySFS=/u/home/a/ab08028/klohmueldata/annabel_data/bin/easySFS/easySFS.abContinueMod.py 

gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/
scriptdir=${gitdir}/analyses/generate_sfs/


easySFS=$scriptdir/easySFS.abModified.2.noInteract.Exclude01Sites.20181121.py  # this is my modification
# this version of script excludes sites that are 0-1 across all populations (maybe) -- not sure if it does yet. 

## choose your projections:
#projections="20,20,20,20,20"
### NOTE: projection values must be in same order as populations are in your popFile (this isn't ideal -- at some point I am going to modify the easySFS script)
# note that order is CA,AK,AL,COM,KUR 

outdir=/u/flashscratch/a/ab08028/captures/analyses/SFS/$genotypeDate/easySFS
mkdir -p $outdir
# had to modify easySFS so that it wouldn't prompt a "yes/no" response about samples that are missing from VCF file

# make sure vcf isn't zipped

snpvcf=neutral_snp_8b_forEasySFS_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf
allSitesvcf=neutral_all_8_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters_raw_variants.vcf.gz
#vcf='neutral.all_8a_rmRelatives_rmAdmixedOutliers_passingBespoke_maxNoCallFrac_'${noCallFrac}'_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf'
# hoping I can get it to work with this neutral all-sites file so that I have monomorphic sites sampled as well. This would be ideal.
#gunzip ${vcf}.gz

### NOTE: projection values must be in same order as populations are in your popFile (this isn't ideal -- at some point I am going to modify the easySFS script)
# note that order is CA,AK,AL,COM,KUR 
$easySFS -i $vcfdir/${snpvcf} -p $popFile -a -v --proj $projections -f -o $outdir/projection-${todaysdate}
# test run only had --proj 20 (maybe just projects 1 population? will have to see)
# -f forces overwrite of outdir
$bgzip ${vcf}
# then do for SYN and MIS (eventually)
########## get counts of monomorphic sites to add to the SFSes ############
python $scriptdir/getMonomorphicProjectionCounts.py --vcf $vcfdir/${allSitesvcf} --popMap $popFile --proj $projections --popIDs CA,AK,AL,COM,KUR --outdir $outdir/projection-${todaysdate}

# need to add these to the 0/0 bin for everything. How to do that? Manually? dadi could do in python when I parse results because it's masked anyway. for fastsimcoal need to add to SFS directly 
# in the 0/0 bin. (but want to add to what's already there.)
# maybe a little R or python script? (coudl read in the fsc SFSes and add this to 0/0 bin?)
########### eventually; to start with can do by hand #######


######## run a little test:
# script=/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/analyses/generate_sfs/getMonomorphicProjectionCounts.py
# python $scriptdir/$script --vcf "/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_filtering/sandbox/dummyVCF.forSandbox.allSites_5_passingFilters.vcf.gz" \
# --popMap /Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/information/samples/easySFSPopMapFiles/samplesPop.Headers.forEasySFS.2.goingtoModifyFor20181119.txt \
# --proj 14,16,16,16,14 --popIDs CA,AK,AL,COM,KUR \
# --outdir /Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/troubleshooting/testRunsEasySFSModificatoins