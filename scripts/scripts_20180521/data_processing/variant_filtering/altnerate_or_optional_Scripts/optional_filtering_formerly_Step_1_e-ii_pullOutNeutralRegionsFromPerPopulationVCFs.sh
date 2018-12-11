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
module load bedtools
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
suffix=all_9_rmAllHet_rmRelativesAdmixed_passingAllFilters_allCalled # vcf suffix; all sites called in all inds. sites that are hets in all inds removed
# 20181018 bed file of neutral sites that have been called (min 10kb from genes, not in CpG island, doesn't blast to zebra fish); max no call frac is 0.9 (liberal)
neutralBed=${vcfdir}/bedCoords/all_8_rmRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_0.9_rmBadIndividuals_passingFilters.min10kb.fromExon.noCpGIsland.noRepeat.noFish.0based.sorted.merged.useThis.bed

populations="CA AK AL COM KUR"
# intersect neutral bed with population vcf files that have all sites called in all individuals
for pop in $populations
do
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/populationVCFs/${pop}_${suffix}.vcf.gz \
-o $outdir/${pop}_neutral_${suffix}.vcf.gz \
-L $neutralBed
done

######### For records, get a bed file of the coords that overlap between pop vcfs and neutBed, and the total amnt of neut sequence per pop: 
echo -e 'pop\ttotalCalledNeutralSites' > ${vcfdir}/bedCoords/neutralCallableSites_perPop/summary.neutralCallableSites.perPop.txt
for pop in $populations
do
echo $pop
# for reference, want to get bed file of callable sites in neutral regions per population
bedtools intersect -a ${vcfdir}/populationVCFs/${pop}_${suffix}.vcf.gz -b $neutralBed | awk '{OFS="\t"; print $1,$2-1,$2}' | sort -k1,1 -k2,2n | bedtools merge -i stdin > ${vcfdir}/bedCoords/neutralCallableSites_perPop/${pop}_${suffix}.neutral.callableSites.0based.bed
# and the total amount of neutral sequence:
totalNeut=`awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' ${vcfdir}/bedCoords/neutralCallableSites_perPop/${pop}_${suffix}.neutral.callableSites.0based.bed`
echo -e ${pop}'\t'${totalNeut} >> ${vcfdir}/bedCoords/neutralCallableSites_perPop/summary.neutralCallableSites.perPop.txt
totalNeut=''
done

