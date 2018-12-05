source /u/local/Modules/default/init/modules.sh
module load java
module load bedtools
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta

snpNoCallFrac=0.2
genotypeDate=20181119

################ Just Baja - California Structure ##################
# select 2 california individuals at random (a few selections)
# not RWAB; not low coverage
# choose from: 
# 114_Elut_CA_214
# 115_Elut_CA_305
# 116_Elut_CA_307
# 117_Elut_CA_315
# 139_Elut_CA_390
# 140_Elut_CA_403
# 141_Elut_CA_419
# 142_Elut_CA_365 
# these aren't low coverage or RWAB -- chose 6 (3 pairs)
ind1=114_Elut_CA_214
ind2=117_Elut_CA_315

ind3=140_Elut_CA_403
ind4=142_Elut_CA_365

ind5=115_Elut_CA_305
ing6=139_Elut_CA_390

# temprorarly using QD2 dir (change once my script finishes)
vcfdir=/u/flashscratch/a/ab08028/captures/vcf_filtering/${genotypeDate}_filtered/QD2Filter_deleteAfterCompare ## temporary!!! 
mkdir -p $vcfdir/subsampledVCFs
vcf=snp_7_maxNoCallFrac_${snpNoCallFrac}_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf.gz
# pair these with the Baja for PCA: 

# first pair: 
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant $vcfdir/$vcf \
-o ${vcfdir}/subsampledVCFs/subsample.1.BAJ.${ind1}.{$ind2}.$vcfFile \
-se $ind1 \
-se $ind2 \
-se '.+_ELUT_BAJ_.+' 

# second pair: 
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant $vcfdir/$vcfFile \
-o ${vcfdir}/subsampledVCFs/subsample.2.BAJ.${ind3}.{$ind4}.$vcfFile \
-se $ind3 \
-se $ind4 \
-se '.+_ELUT_BAJ_.+' 

# third pair:
# second pair: 
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant $vcfdir/$vcfFile \
-o ${vcfdir}/subsampledVCFs/subsample.3.BAJ.${ind3}.{$ind4}.$vcfFile \
-se $ind5 \
-se $ind6 \
-se '.+_ELUT_BAJ_.+' 