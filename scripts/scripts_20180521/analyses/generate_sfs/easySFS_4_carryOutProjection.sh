#! /bin/bash
#$ -cwd
#$ -l h_rt=10:00:00,h_data=20G,highp
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
module load R

bgzip=/u/home/a/ab08028/klohmueldata/annabel_data/bin/tabix-0.2.6/bgzip
todaysdate=`date +%Y%m%d`
genotypeDate=20181119
noCallFrac=1.0
vcfdir=/u/flashscratch/a/ab08028/captures/vcf_filtering/${genotypeDate}_filtered/neutral_and_cds_VCFs/
popFile=/u/flashscratch/a/ab08028/captures/samples/samplesPop.Headers.forEasySFS.3.20181119.txt
allSamplesHetFilter=0.75 # het filtering done across all samples
perPopHetFilter=0.75 # the vcf file has already had some degree of maxHetFiltering. now per-population, easy sfs will do per-population filtering at this level  


gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/
scriptdir=${gitdir}/analyses/generate_sfs/


easySFS=$scriptdir/easySFS.abModified.3.noInteract.Exclude01Sites.HetFiltering.20181121.py  # this is my modification
# this version of script filters het sites that are excessively heterozygous 

## choose your projections: choosing for now: 
# CA,AK,AL,COM,KUR 
# choosing projections to maximize snps:
#Choices based on 20181220 projection preview: 
#AK Ð 14 (max)
#AL Ð 20 (max)
#CA Ð 12 (max would be at 8)
#COM Ð 34 (max)
#KUR Ð 12 (max would be at 10)
 #CA,AK,AL,COM,KUR 
projections="12,14,20,34,12" # updated these values on 20181220
### NOTE: projection values must be in same order as populations are in your popFile (this isn't ideal -- at some point I am going to modify the easySFS script)
# note that order is CA,AK,AL,COM,KUR 

############################### NEUTRAL SITE PROJECTIONS #############################

outdir=/u/flashscratch/a/ab08028/captures/analyses/SFS/$genotypeDate/easySFS/neutral/projection-${todaysdate}
mkdir -p $outdir
# had to modify easySFS so that it wouldn't prompt a "yes/no" response about samples that are missing from VCF file
# write projection choices into a readme

echo "CA,AK,AL,COM,KUR : $projections " > $outdir/projectionChoices.${todaysdate}.txt
# make sure vcf isn't zipped

allVCF=neutral.all_9_maxHetFilter_${allSamplesHetFilter}_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters_raw_variants.vcf.gz
snpVCF=neutral.snp_9b_forEasySFS_maxHetFilter_${allSamplesHetFilter}_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf


### NOTE: projection values must be in same order as populations are in your popFile (this isn't ideal -- at some point I am going to modify the easySFS script)
# note that order is CA,AK,AL,COM,KUR 
$easySFS -i $vcfdir/neutralVCFs/${snpVCF} -p $popFile -a -v --proj $projections -f -o $outdir -maxHetFilter $perPopHetFilter
# test run only had --proj 20 (maybe just projects 1 population? will have to see)
# -f forces overwrite of outdir
# $bgzip ${vcf}
# then do for SYN and MIS (eventually)
########## get counts of monomorphic sites to add to the SFSes ############
python $scriptdir/getMonomorphicProjectionCounts.py --vcf $vcfdir/neutralVCFs/${allVCF} --popMap $popFile --proj $projections --popIDs CA,AK,AL,COM,KUR --outdir $outdir

###### you are here in interactive node ########
############## adding monomorphic sites to fsc SFSes #####################
 # this script add monomorphic sites to 0 bin of fsc sfses. doesn't add them to dadi SFSes because those sites are masked anyway. 
Rscript $scriptdir/easySFS_5_addInMonomorphicSites.R --dataDir $outdir --popFile $popFile # will write them out in your data dir in new directories





##################################### coding sites #############################
cdsVCF=cds_${allVCF}
synVCF=syn_vep_cds_${snpVCF%.gz}
misVCF=missense_vep_cds_${snpVCF%.gz}

outdir=/u/flashscratch/a/ab08028/captures/analyses/SFS/$genotypeDate/easySFS/cds/synonymous/projection-${todaysdate}
mkdir -p $outdir
$easySFS -i $vcfdir/cdsVCFs/${synVCF} -p $popFile -a -v --proj $projections -f -o $outdir -maxHetFilter $perPopHetFilter

outdir=/u/flashscratch/a/ab08028/captures/analyses/SFS/$genotypeDate/easySFS/cds/missense/projection-${todaysdate}
mkdir -p $outdir
$easySFS -i $vcfdir/cdsVCFs/${misVCF} -p $popFile -a -v --proj $projections -f -o $outdir -maxHetFilter $perPopHetFilter


########## get counts of monomorphic cds sites (but don't add to SFSes yet; need to think about scaling more)############
outdir=/u/flashscratch/a/ab08028/captures/analyses/SFS/$genotypeDate/easySFS/cds/monomorphicCDSSites/projection-${todaysdate}
mkdir -p $outdir
python $scriptdir/getMonomorphicProjectionCounts.py --vcf $vcfdir/cdsVCFs/${cdsVCF} --popMap $popFile --proj $projections --popIDs CA,AK,AL,COM,KUR --outdir $outdir


############ troubleshooting het filter: 
##outdir=/u/flashscratch/a/ab08028/captures/analyses/SFS/$genotypeDate/easySFS/troubleshoot-hetFilter
##mkdir -p $outdir
#maxHetFilter=0.75 # this is now per-population level filtering of hets.
# $scriptdir/easySFS.abModified.3.noInteract.Exclude01Sites.HetFiltering.20181121.py -i $vcfdir/${snpvcf} -p $popFile -a -v --proj $projections -f -o $outdir -maxHetFilter $maxHetFilter


# get DP dist of sites that pass/fail a het filter:
#wd=/u/flashscratch/a/ab08028/sandbox/hetFilters
#script=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/analyses/generate_sfs/sandbox.getDP.QD.QUAL.easySFS.abModified.2.noInteract.Exclude01Sites.HetFilterExperiments.20181121.py
#maxHetFilter=0.8
#python $script -i $vcfdir/${snpvcf} -p $popFile -a -v --proj $projections -f -o $outdir -maxHetFilter $maxHetFilter -dpFile $wd/DP.QD.QUAL.dist.HetFilter.${maxHetFilter}.txt