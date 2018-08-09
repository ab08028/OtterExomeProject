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
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar

outdir=$wd/${rundate}_filtered # date you called genotypes
#################################################################################
############################ Separate by population #############################
#################################################################################
# want to separate vcfs by population
mkdir $outdir/populationVCFs
populations="CA AL KUR AK"
# treat commanders specially

# need to combine Bering and Medny separately (and label as COM (commanders))
########### bering and medny: ###########
# note: I'm making it so that every individual in the population has to be called at the site,
# but across all the populations only 80% have to be called at the site. I think this makes sense

# do snps: [fast ~5-10 min]
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${outdir}/'snp_7_maxNoCallFrac_'${noCallFrac}'_passingBespoke_passingAllFilters_postMerge_'${infile} \
-o ${outdir}/populationVCFs/COM_'snp_7_allIndCalled_passingBespoke_passingAllFilters_postMerge_'${infile} \
-se '.+_Elut_BER_.+' \
-se '.+_Elut_MED_.+' \
--maxNOCALLfraction 0
# project down SFS? or not? <--- question to ask people. 

# do nv:
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${outdir}/'nv_7_maxNoCallFrac_'${noCallFrac}'_passingBespoke_passingAllFilters_postMerge_'${infile} \
-o ${outdir}/populationVCFs/COM_'nv_7_allIndCalled_passingBespoke_passingAllFilters_postMerge_'${infile} \
-se '.+_Elut_BER_.+' \
-se '.+_Elut_MED_.+' \
--maxNOCALLfraction 0

# go through rest of populations: (this assumes you've already removed any bad individuals during filtering)
for pop in $populations
do
echo $pop

# do snps: [fast ~5-10 min]
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${outdir}/'snp_7_maxNoCallFrac_'${noCallFrac}'_passingBespoke_passingAllFilters_postMerge_'${infile} \
-o ${outdir}/populationVCFs/${pop}_'snp_7_maxNoCallFrac_'${noCallFrac}'_passingBespoke_passingAllFilters_postMerge_'${infile} \
-se '.+_Elut_${pop}_.+'


# do nv:
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${outdir}/'nv_7_maxNoCallFrac_'${noCallFrac}'_passingBespoke_passingAllFilters_postMerge_'${infile} \
-o ${outdir}/populationVCFs/${pop}_'nv_7_maxNoCallFrac_'${noCallFrac}'_passingBespoke_passingAllFilters_postMerge_'${infile} \
-se '.+_Elut_${pop}_.+'

done
########### test:

# does this update AN fields? check:
pop=alaska
popHeaderDir=/u/flashscratch/a/ab08028/captures/samples/popHeaderFiles/popHeaders_allInd/
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant /u/flashscratch/a/ab08028/captures/vcf_filtering/20180724_filtered/snp_5_passingAllFilters_postMerge_raw_variants.vcf.gz \
-o /u/flashscratch/a/ab08028/captures/vcf_filtering/20180724_filtered/CA_snp_5_passingAllFilters_postMerge_raw_variants.vcf.gz \
-sn *Elut_CA_* \
--maxNOCALLfraction 0
# currently doesn't output any snps...?