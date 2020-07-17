#! /bin/bash
#$ -cwd
#$ -l h_rt=200:00:00,h_data=28G,highp,arch=intel*
#$ -N jointGenoSSO
#$ -o /u/flashscratch/a/ab08028/captures/reports/GATK/
#$ -e /u/flashscratch/a/ab08028/captures/reports/GATK/
#$ -m abe
#$ -M ab08028
 
rundate=`date +%Y%m%d`  # date you call genotypes; can use it to distinguish from multiple times you call gts

source /u/local/Modules/default/init/modules.sh
module load java/1.8.0_111
module load samtools

# file locations:
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures
# note: covered intervals already taken into account in hapcaller,so don't need to use -L flag here
# ferret reference:
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/sea_otter_genome/dedup_99_indexed_USETHIS/sea_otter_23May2016_bS9RH.deduped.99.fasta
REFPREFIX=sea_otter_23May2016_bS9RH.deduped.99
REFSHORTCODE=SSO

GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar
outdir=$wd/vcfs/vcf_${rundate}_${REFSHORTCODE}

mkdir -p $outdir
bgzip=~/klohmueldata/annabel_data/bin/tabix-0.2.6/bgzip

# get vcfRunningList:
# manually exclude any low coverage samples you want 
# this will joint call everything that's in the gVCFs directory (careful with that)
> $wd/samples/vcfRunningList.${rundate}.${REFSHORTCODE}.list
samples=`ls $wd/gvcfs/$REFSHORTCODE/*gz`
for i in $samples
do
echo $i >> $wd/samples/vcfRunningList.${rundate}.${REFSHORTCODE}.list # this list will record everything that went into genotyping on XX date.
done
# This is for all populations
# for now not changing het, or other parameters. 
# vcf list must not have any comment lines; must end in .list; can have full paths
java -jar $GATK \
  -T GenotypeGVCFs \
  -R $REFERENCE \
  -V $wd/samples/vcfRunningList.${rundate}.${REFSHORTCODE}.list \
  --includeNonVariantSites \
  -o ${outdir}/raw_variants.vcf.gz
 
# zip the file if GATK didn't (may not be necessary)
### $bgzip ${OUT_DIR}/elut.raw_variants.${rundate}.vcf
# for some reason an index isn't always made? 
  
# tabix -p vcf ${OUT_DIR}/elut.raw_variants.${rundate}.vcf.gz
