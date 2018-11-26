#! /bin/bash
#$ -cwd
#$ -l h_rt=10:00:00,h_data=16G,highp
#$ -N easySFSFilePrep
#$ -o /u/flashscratch/a/ab08028/captures/reports/SFS
#$ -e /u/flashscratch/a/ab08028/captures/reports/SFS
#$ -m abe
#$ -M ab08028
########## Prepare files for EasySFS
### Need: neutral regions (all and snps only) (leave VEP for later)

source /u/local/Modules/default/init/modules.sh
module load java
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar

#### parameters:
rundate=20180806 # date genotypes were called (vcf_20180806 includes capture 02)
noCallFrac=0.9 # maximum fraction of genotypes that can be "no call" (./.) that was used in previous steps in previous "all sites" files (lenient cutoff)
#snpNoCallFrac=0.2 # max missingness allowed in snp file (stricter cutoff)
perPopNoCallFrac=0 # max missingness allowed in final file for each pop for sfs (super strict cutoff)

#### file locations
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/vcf_filtering
mkdir -p $wd
infile=raw_variants.vcf.gz ### make sure this doesn't have a path as part of its name! just infile names
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta
vcfdir=$wd/${rundate}_filtered # date you called genotypes
neutBed=${vcfdir}/bedCoords/all_8_rmRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_0.9_rmBadIndividuals_passingFilters.min10kb.fromExon.noCpGIsland.noRepeat.noFish.0based.sorted.merged.useThis.bed

mkdir -p $vcfdir/populationVCFs
mkdir -p $vcfdir/populationVCFs/admixedVCFs

# this vcf is the result of STep 1 d 
# eventually want to use all_8a ...rmAdmixedOutliers ... 
# so that I have monomorphic info
# but for now going to use the snp file 
vcf='all_8a_rmRelatives_rmAdmixedOutliers_passingBespoke_maxNoCallFrac_'${noCallFrac}'_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf.gz'

#vcf='snp_8a_rmRelatives_rmAdmixedOutliers_passingBespoke_maxNoCallFrac_'${noCallFrac}'_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf.gz'
mkdir -p $vcfdir/neutralVCF_allPops

java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/$vcf \
-o $vcfdir/neutralVCF_allPops/neutral.${vcf} \
-L $neutBed
# don't output as bgzipped 

# then do this for VEP output.
# I think I want to run VEP on all-sites and then sep by syn/mis etc. <-- this needs some figuring out.

#### eventually going to do this for all_8 -- but that will take time to prepare.