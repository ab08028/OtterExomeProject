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

snpvcf=neutral.snp_8b_forEasySFS_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf
allSitesvcf=neutral.all_8_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters_raw_variants.vcf.gz
#vcf='neutral.all_8a_rmRelatives_rmAdmixedOutliers_passingBespoke_maxNoCallFrac_'${noCallFrac}'_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf'
# hoping I can get it to work with this neutral all-sites file so that I have monomorphic sites sampled as well. This would be ideal.
#gunzip ${vcf}.gz

### NOTE: projection values must be in same order as populations are in your popFile (this isn't ideal -- at some point I am going to modify the easySFS script)
# note that order is CA,AK,AL,COM,KUR 
$easySFS -i $vcfdir/${snpvcf} -p $popFile -a -v --proj $projections -f -o $outdir
# test run only had --proj 20 (maybe just projects 1 population? will have to see)
# -f forces overwrite of outdir
# $bgzip ${vcf}
# then do for SYN and MIS (eventually)
########## get counts of monomorphic sites to add to the SFSes ############
python $scriptdir/getMonomorphicProjectionCounts.py --vcf $vcfdir/${allSitesvcf} --popMap $popFile --proj $projections --popIDs CA,AK,AL,COM,KUR --outdir $outdir
# test: started 12:20
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

## testing that my monormophic script matches what GATK would give:
# zcat $vcfdir/$allSitesvcf | grep -v "#" | grep "0/0" | grep -v "0/1" | grep -v "1/1" -c
# # lines that are all 0/0 and no 0/1 and 1/1 (doesn't account for missing data or projection values yet: 6585427
# # Okay, so they are on the right order of magnitude. 
# # how to test that this is right? compare to GATK with nocallfrac and restrict to just one population
# vcfdir=/u/flashscratch/a/ab08028/captures/vcf_filtering/20181119_filtered/neutral_and_cds_VCFs/neutralVCFs/
# allSitesvcf=neutral.all_8_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters_raw_variants.vcf.gz
# wd=/u/flashscratch/a/ab08028/sandbox/monomorphicCounts
# GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar
# REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta
# java -jar $GATK \
# -R $REFERENCE \
# -T SelectVariants \
# --variant ${vcfdir}/$allSitesvcf \
# --maxNOCALLnumber 6 \
# -selectType NO_VARIATION \
# -se 'Elut_CA' \
# -o $wd/test.MonomorphicCounter.CA.vcf.gz
# runs in 15 min
# # line count (without # lines) should equal 5790572 for CA count.
outdir=/u/flashscratch/a/ab08028/captures/analyses/SFS/$genotypeDate/easySFS/troubleshoot-hetFilter
$scriptdir/sandbox.easySFS.abModified.2.noInteract.Exclude01Sites.HetFilterExperiments.20181121.py -i $vcfdir/${snpvcf} -p $popFile -a -v --proj $projections -f -o $outdir
