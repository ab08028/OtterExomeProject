#! /bin/bash
#$ -cwd
#$ -l h_rt=50:00:00,h_data=16G,highp,arch=intel*
#$ -N 1b_vcf_filtering
#$ -o /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -e /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -m abe
#$ -M ab08028

rundate=20180806
#### file locations
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/vcf_filtering
infile=raw_variants.vcf.gz ### make sure this doesn't have a path as part of its name! just infile names
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta

# need to remove bad individuals from lists manually before this 
popHeaderDir=$SCRATCH/captures/samples/popHeaders_rmBadInd/

outdir=$wd/${rundate}_filtered # date you called genotypes
#################################################################################
############################ Separate by population #############################
#################################################################################
# want to separate vcfs by population
mkdir $outdir/populationVCFs
populations="alaska aleutian commanders california kuril"
for pop in $populations
do
echo $pop

# do snps:
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${outdir}/'snp_7_80perc_passingBespoke_passingAllFilters_postMerge_'${infile} \
-o ${outdir}/populationVCFs/${pop}_'snp_7_80perc_passingBespoke_passingAllFilters_postMerge_'${infile} \
-IDs ${popHeaderDir}/${pop}.txt

# do nv:
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${outdir}/'nv_7_80perc_passingBespoke_passingAllFilters_postMerge_'${infile} \
-o ${outdir}/populationVCFs/${pop}_'nv_7_80perc_passingBespoke_passingAllFilters_postMerge_'${infile} \
-IDs ${popHeaderDir}/${pop}.txt

done

# does this update AN fields? check:
pop=california

java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant /u/flashscratch/a/ab08028/captures/vcf_filtering/20180724_filtered/snp_5_passingAllFilters_postMerge_raw_variants.vcf.gz
-o /u/flashscratch/a/ab08028/captures/vcf_filtering/20180724_filtered/${pop}_snp_5_passingAllFilters_postMerge_raw_variants.vcf.gz \
-IDs ${popHeaderDir}/${pop}.txt