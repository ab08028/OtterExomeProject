#! /bin/bash
#$ -cwd
#$ -l h_rt=50:00:00,h_data=16G,highp,arch=intel*
#$ -N vcf1e_filtering
#$ -o /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -e /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -m abe
#$ -M ab08028

# modules
source /u/local/Modules/default/init/modules.sh
module load java
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar

rundate=20180806
noCallFrac=0.2
#### file locations
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/vcf_filtering
infile=raw_variants.vcf.gz ### make sure this doesn't have a path as part of its name! just infile names
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta

vcfdir=$wd/${rundate}_filtered # date you called genotypes

# neutral sites that have been called (min 10kb from genes, not in CpG island, doesn't blast to zebra fish)
neutralBed=${vcfdir}/bedCoords/all_7_passingBespoke.min10kb.fromExon.noCpGIsland.noRepeat.noFish.0based.sorted.merged.useThis.bed


populations="CA AL KUR AK COM"

for pop in $populations
do
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/populationVCFs/${pop}_'all_7_passingAllFilters_allCalled.vcf.gz' \
-o ${vcfdir}/populationVCFs/${pop}_'neutral_7_passingAllFilters_allCalled.vcf.gz' \
-L $neutralBed
done
