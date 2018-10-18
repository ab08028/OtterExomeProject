#! /bin/bash
#$ -cwd
#$ -l h_rt=50:00:00,h_data=16G,highp,arch=intel*
#$ -N vcf1e_filtering
#$ -o /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -e /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -m abe
#$ -M ab08028

# 20181018 decided this is an important step (Tanya's scripts for pi/total snps dont use chr designations so don't work for my data; remaking my scripts)
# modules
source /u/local/Modules/default/init/modules.sh
module load java
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar

rundate=20180806
noCallFrac=0.2
#### file locations
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/vcf_filtering
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta

vcfdir=$wd/${rundate}_filtered # date you called genotypes
outdir=$vcfdir/populationVCFs/neutralVCFs
mkdir -p $outdir
suffix=all_9_rmAllHet_rmRelativesAdmixed_passingAllFilters_allCalled # vcf suffix
# 20181018 neutral sites that have been called (min 10kb from genes, not in CpG island, doesn't blast to zebra fish)
neutralBed=${vcfdir}/bedCoords/all_8_rmRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_0.9_rmBadIndividuals_passingFilters.min10kb.fromExon.noCpGIsland.noRepeat.noFish.0based.sorted.merged.useThis.bed

populations="CA AK AL COM KUR"

for pop in $populations
do
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/populationVCFs/${pop}_${suffix}.vcf.gz \
-o $outdir/${pop}_neutral_${suffix}.vcf.gz \
-L $neutralBed
done
