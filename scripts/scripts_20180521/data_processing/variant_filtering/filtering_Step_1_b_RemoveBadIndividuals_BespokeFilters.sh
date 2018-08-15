#! /bin/bash
#$ -cwd
#$ -l h_rt=50:00:00,h_data=16G,highp,arch=intel*
#$ -N vcf1b_filtering
#$ -o /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -e /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -m abe
#$ -M ab08028

# modules
source /u/local/Modules/default/init/modules.sh
module load java
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar
tabix=/u/home/a/ab08028/klohmueldata/annabel_data/bin/tabix-0.2.6/tabix

#### parameters:
rundate=20180806 # date genotypes were called (vcf_20180806 includes capture 02)
noCallFrac=0.2 # maximum fraction of genotypes that can be "no call" (./.) # note that genotypes can still "PASS" if 


#### file locations
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/vcf_filtering
mkdir -p $wd
indir=$SCRATCH/captures/vcfs/vcf_${rundate}
infile=raw_variants.vcf.gz ### make sure this doesn't have a path as part of its name! just infile names
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta
scriptdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_filtering

bespokeFilterScript=$scriptdir/filtering_bespokeFiltersAndChecks.py

badIndFile=$SCRATCH/samples/bad.individuals.toRemove.DuringFiltering.txt # IDs of individuals to remove (get from step 1a)

# location of vcf checking and filtering script
# incompatible scaffolds: repeatMaskCoords=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/repeatMaskingCoordinates/masking_coordinates.bed
# repeat masking is optional: my target captrue was designed away from repeat regions, so not a huge deal 
# if you don't include this; but trying to be extra thorough.
# these coords are downloaded from NCBI FTP site for mustela putorius furo ; if using elut genome, use Annotation repeat coords; 
# ftp://ftp.ncbi.nlm.nih.gov/genomes/Mustela_putorius_furo/ <-- location of mustela files
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Mustela_putorius_furo/masking_coordinates.gz

outdir=$wd/${rundate}_filtered # date you called genotypes
mkdir -p $outdir


######################## Remove bad individuals (relatives, low coverage, PCA outliers, admixed, etc.) ###############
# which means that you'll rerun Step1b multiple times most likely (but hopefully DON'T have to rerun step 1a)
# bad individuals: add/remove as needed based on how many missing genotypes there are per individual (>15 million, cut)
ind1=
ind2=
ind3=
ind4=
ind5=
ind6=
ind7=
ind8=
ind9=
ind10=
ind11=

echo "starting step 6: remove bad individuals"
echo "These are the bad (high dup %, high missing GTs, low coverage) individuals that are getting removed: $ind1 $ind2 $ind3 $ind4 $ind5 $ind6 $ind7 $ind8 $ind9 $ind10 $ind11"
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${outdir}/'all_5_passingFilters_'${infile} \
-o ${outdir}/'all_6_rmBadIndividuals_passingFilters_'${infile} \
-xl_sn ${ind1} \
-xl_sn $ind2 \
-xl_sn $ind3 \ 
-xl_sn $ind4 \
-xl_sn $ind5 \
-xl_sn $ind6 \
-xl_sn $ind7 \
-xl_sn $ind8 \ 
-xl_sn $ind9 \
-xl_sn ${ind10} \
-xl_sn ${ind11}

 
# note that 4 extremely low samples were already excluded prior to genotype calling from the medgenome captures
# note I need to add the other bad individuals from 
echo "done with step 6: remove bad individuals"
#################################################################################
############################ RUN BESPOKE FILTERS and UPDATE AN/AC ##########################
#################################################################################
# These filters will check for:
# usage: python [script] [full path to invcf] [full path to out vcf] [full path to error file] [max no call fraction]
# 1. ref or alt alleles that aren't a single letter (AGCT) or . (alt)
# 2. genotypes that aren't in 0/0, 0/1, 1/1 or ./. (maybe it's phased, etc)
# 3. must have qual score
# 4. must be PASS for site
# 5. make sure not missing DP, AN, GT, AD, DP or GQ/RGQ
# 6. make sure no called genotype is missing any info from the genotype info field
# 7. gets rid of sites where all calls are 0/1 (all hets)
# 8. updates AN and AC based on final sets of calls (these aren't updated when GATK does genotype filtering)
# 9. filters sites that exceed the maximum no-call fraction (supplied by user)
# this script does NOT: change any genotypes; do any genotype filtering; change any FT fields for genotypes (./. gts will still be PASS if they started as ./. -- bit of GATK weirdness that isn't fatal)
echo "starting step 7a: carrying out bespoke filtering (includes filtering out sites exceeding max no-call frac: $noCallFrac)"

python $bespokeFilterScript ${outdir}/'all_6_rmBadIndividuals_passingFilters_'${infile} \
${outdir}/'all_7_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile%.gz} \
${outdir}/'fail_all_7_FAILINGBespoke_passingFilters_'${infile%.vcf.gz}.txt \
$noCallFrac
# gzip the result:
gzip  ${outdir}/'all_7_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile%.gz}
$tabix -p ${outdir}/'all_7_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} # index the vcf

echo "done with step 7a: carrying out bespoke filtering"

#################################################################################
############################ Get final sets of variant / invariant###############
#################################################################################
echo "starting step 7b: select final passing snps from merged file"
# Some sites that may have started as variant may have become INVARIANT by the end.
# Want these to end up in the INVARIANT category. 
# select the biallelic snps: 
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-V ${outdir}/'all_7_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile%.gz} \
--restrictAllelesTo BIALLELIC \
--selectTypeToInclude SNP \
-o ${outdir}/'snp_7_maxNoCallFrac_'${noCallFrac}'_passingBespoke_passingAllFilters_postMerge_'${infile}
echo "done step 7b: select final passing snps from merged file"

## Select the invariants:
echo "starting step 6c: select final passing nonvariant sites from merged file"

java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-V ${outdir}/'all_7_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile%.gz} \
--selectTypeToInclude NO_VARIATION \
--selectTypeToExclude INDEL \
--maxNOCALLfraction $noCallFrac \
-o ${outdir}/'nv_7_maxNoCallFrac_'${noCallFrac}'_passingBespoke_passingAllFilters_postMerge_'${infile}
echo "done step 7c: select final passing nonvariant sites from merged file"


