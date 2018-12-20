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
#todaysdate=20181212
todaysdate=`date +%Y%m%d`
genotypeDate=20181119
noCallFrac=1.0
vcfdir=/u/flashscratch/a/ab08028/captures/vcf_filtering/${genotypeDate}_filtered/neutral_and_cds_VCFs/neutralVCFs/
popFile=/u/flashscratch/a/ab08028/captures/samples/samplesPop.Headers.forEasySFS.3.20181119.txt
# this has admixed in it , but they aren't in pop file
#easySFS=/u/home/a/ab08028/klohmueldata/annabel_data/bin/easySFS/easySFS.abContinueMod.py 

gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/
scriptdir=${gitdir}/analyses/generate_sfs/


easySFS=$scriptdir/easySFS.abModified.2.noInteract.Exclude01Sites.20181121.py  # this is my modification
# this version of script excludes sites that are 0-1 across all populations (maybe) -- not sure if it does yet. 

## choose your projections: choosing for now: 
# CA,AK,AL,COM,KUR 
projections="12,14,16,16,12" # may change these after lab meeting
### NOTE: projection values must be in same order as populations are in your popFile (this isn't ideal -- at some point I am going to modify the easySFS script)
# note that order is CA,AK,AL,COM,KUR 

outdir=/u/flashscratch/a/ab08028/captures/analyses/SFS/$genotypeDate/easySFS/projection-${todaysdate}
mkdir -p $outdir
# had to modify easySFS so that it wouldn't prompt a "yes/no" response about samples that are missing from VCF file
# write projection choices into a readme

echo "CA,AK,AL,COM,KUR : $projections " > $outdir/projectionChoices.${todaysdate}.txt
# make sure vcf isn't zipped

allVCF=all_9_maxHetFilter_${maxHetFilter}_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters_raw_variants.vcf.gz
snpVCF=snp_9b_forEasySFS_maxHetFilter_${maxHetFilter}_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf.gz


### NOTE: projection values must be in same order as populations are in your popFile (this isn't ideal -- at some point I am going to modify the easySFS script)
# note that order is CA,AK,AL,COM,KUR 
$easySFS -i $vcfdir/${snpvcf} -p $popFile -a -v --proj $projections -f -o $outdir
# test run only had --proj 20 (maybe just projects 1 population? will have to see)
# -f forces overwrite of outdir
# $bgzip ${vcf}
# then do for SYN and MIS (eventually)
########## get counts of monomorphic sites to add to the SFSes ############
python $scriptdir/getMonomorphicProjectionCounts.py --vcf $vcfdir/${allSitesvcf} --popMap $popFile --proj $projections --popIDs CA,AK,AL,COM,KUR --outdir $outdir

############ troubleshooting: 
outdir=/u/flashscratch/a/ab08028/captures/analyses/SFS/$genotypeDate/easySFS/troubleshoot-hetFilter
mkdir -p $outdir
maxHetFilter=0.5
$scriptdir/sandbox.easySFS.abModified.2.noInteract.Exclude01Sites.HetFilterExperiments.20181121.py -i $vcfdir/${snpvcf} -p $popFile -a -v --proj $projections -f -o $outdir -maxHetFilter $maxHetFilter


# get DP dist of sites that pass/fail a het filter:
wd=/u/flashscratch/a/ab08028/sandbox/hetFilters
script=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/analyses/generate_sfs/sandbox.getDP.QD.QUAL.easySFS.abModified.2.noInteract.Exclude01Sites.HetFilterExperiments.20181121.py
maxHetFilter=0.8
python $script -i $vcfdir/${snpvcf} -p $popFile -a -v --proj $projections -f -o $outdir -maxHetFilter $maxHetFilter -dpFile $wd/DP.QD.QUAL.dist.HetFilter.${maxHetFilter}.txt