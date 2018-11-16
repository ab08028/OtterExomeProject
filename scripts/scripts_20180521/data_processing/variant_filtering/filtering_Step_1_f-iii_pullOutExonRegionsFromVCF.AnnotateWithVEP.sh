#! /bin/bash
#$ -cwd
#$ -l h_rt=50:00:00,h_data=16G,highp
#$ -N exonVEP
#$ -o /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -e /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -m abe
#$ -M ab08028
#$ -pe shared 3

# 20181018 decided this is an important step (Tanya's scripts for pi/total snps dont use chr designations so don't work for my data; remaking my scripts)
# modules
source /u/local/Modules/default/init/modules.sh
module load java
module load bedtools
module load perl/5.10.1
module load htslib

GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar
vepdir=/u/home/a/ab08028/klohmueldata/annabel_data/bin/ensembl-vep/

rundate=20180806
noCallFrac=0.2
#### file locations
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/vcf_filtering
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta
cdsBed=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/MusPutFuro1.0.91.cdsOnly.0based.sorted.merged.bed
# based on MusPutFuro1.0.91.cdsOnly.gff
vcfdir=$wd/${rundate}_filtered # date you called genotypes
outdir=$vcfdir/populationVCFs/cdsVCFs
mkdir -p $outdir
mkdir -p ${vcfdir}/bedCoords/cdsCallableSites_perPop/
suffix=all_9_rmAllHet_rmRelativesAdmixed_passingAllFilters_allCalled # vcf suffix; all sites called in all inds. sites that are hets in all inds removed

populations="CA AK AL COM KUR"
# intersect cds bed with population vcf files that have all sites called in all individuals
for pop in $populations
do
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/populationVCFs/${pop}_${suffix}.vcf.gz \
-o $outdir/${pop}_cds_${suffix}.vcf \
-L $cdsBed
# vcf shouldn't be zipped
# run VEP :
$vepdir/vep -v -i $outdir/${pop}_cds_${suffix}.vcf  \
--cache --force_overwrite --species mustela_putorius_furo \
--variant_class --vcf --canonical \
-o $outdir/${pop}_VEP_cds_${suffix}.vcf


# excluding: numbers, domains
# --fork 3
# select NS (missense) and S in canonical transcripts:
###### S (synonymous):

$vepdir/filter_vep --filter "Consequence is synonymous_variant and CANONICAL is YES" \
--input_file $outdir/${pop}_VEP_cds_${suffix}.vcf \
--output_file $outdir/${pop}_VEP_syn_${suffix}.vcf  \
--force_overwrite


$vepdir/filter_vep --filter "Consequence is missense_variant and CANONICAL is YES" \
--input_file $outdir/${pop}_VEP_cds_${suffix}.vcf \
--output_file $outdir/${pop}_VEP_missense_${suffix}.vcf  \
--force_overwrite

done


## bgzip vcfs at some point
######### For records, get a bed file of the coords that overlap between pop vcfs and neutBed, and the total amnt of neut sequence per pop: 
echo -e 'pop\ttotalCalledcdsSites' > ${vcfdir}/bedCoords/cdsCallableSites_perPop/summary.cdsCallableSites.perPop.txt
for pop in $populations
do
echo $pop
for file in $outdir/${pop}_VEP_cds_${suffix}.vcf $outdir/${pop}_VEP_syn_${suffix}.vcf $outdir/${pop}_VEP_missense_${suffix}.vcf
do
# for reference, want to get bed file of callable sites in cds regions per population
bedtools intersect -a ${vcfdir}/populationVCFs/${pop}_${suffix}.vcf.gz -b $cdsBed | awk '{OFS="\t"; print $1,$2-1,$2}' | sort -k1,1 -k2,2n | bedtools merge -i stdin > ${vcfdir}/bedCoords/cdsCallableSites_perPop/${pop}_${suffix}.cds.callableSites.0based.bed
# and the total amount of coding sequence (cds) sequence:
totalSeq=`awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' ${vcfdir}/bedCoords/cdsCallableSites_perPop/${pop}_${suffix}.cds.callableSites.0based.bed`
echo -e ${pop}'\t'${totalNeut} >> ${vcfdir}/bedCoords/cdsCallableSites_perPop/summary.cdsCallableSites.perPop.txt
totalSeq=''
# then get total Syn and Mis
done
done

# For CA, it's only 12 Mb called. WTF. Where are all my regions?? what is going on?