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
maxHetFilter=0.75 # het filter used across all samples (per population het filter occurs during easy sfs)

genotypeDate=20181119
vcfdir=/u/flashscratch/a/ab08028/captures/vcf_filtering/${genotypeDate}_filtered/neutral_and_cds_VCFs/neutralVCFs
popFile=/u/flashscratch/a/ab08028/captures/samples/samplesPop.Headers.forEasySFS.3.20181119.txt # this doesn't have baja on it; doesn't have any admixed/bad inds on it. 
# this has admixed in it , but they aren't in pop file
easySFS=/u/home/a/ab08028/klohmueldata/annabel_data/bin/easySFS/easySFS.abContinueMod.py
outdir=/u/flashscratch/a/ab08028/captures/analyses/SFS/$genotypeDate/easySFS/projection_preview
mkdir -p $outdir
# had to modify easySFS so that it wouldn't prompt a "yes/no" response about samples that are missing from VCF file
# want it to just output that info and continue , not prompt yes/no.
# this vcf has all snps across all categories (cds, neutral, etc.) with 0.9 max no call frac (v. liberal)
# and has had all individuals removed that won't go into the SFS
# going to do the actual projection for each category of site
#vcf=neutral.snp_8b_forEasySFS_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf
vcf=neutral.snp_9b_maxHetFilter_${maxHetFilter}_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters_raw_variants.vcf

#gunzip $vcfdir/populationVCFs/$vcf it must be unzipped
#( you are here )
$easySFS -i $vcfdir/${vcf} -p $popFile --preview -a -v > $outdir/neutral.snp_9b.easySFS.projPreview.txt

### now plot projections in R and decide on your levels. Actually DO the projections on 

# need to add -a otherwise it selects one snp per locus (was built for RAD data)
# eventually make popFile list consistent with final dataset individuals
# otherwise have to say "yes" on screen
    #Running preview mode. We will print out the results for # of segregating sites
    #for multiple values of projecting down for each population. The dadi
    #manual recommends maximizing the # of seg sites for projections, but also
   	#a balance must be struck between # of seg sites and sample size.

    #For each population you should choose the value of the projection that looks
   # best and then rerun easySFS with the `--proj` flag.
