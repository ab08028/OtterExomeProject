#! /bin/bash
#$ -cwd
#$ -l h_rt=5:00:00,h_data=16G
#$ -N testingGATK
#$ -o /u/flashscratch/a/ab08028/sandbox/reports
#$ -e /u/flashscratch/a/ab08028/sandbox/reports
#$ -m abe
#$ -M ab08028

########## testing filtering with elut reference for now
############ test parameters: 1 chromosome; capture 02 data. #############
module load java
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar

numInd=26
perc80=21
maxNoCall=5
wd=$SCRATCH/sandbox
infile=$wd/elut.raw_variants.20170914.vcf.gz
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/sea_otter_genome/dedup_99_indexed_USETHIS/sea_otter_23May2016_bS9RH.deduped.99.fasta
mychrom=ScbS9RH_6353 # for now  select 1 scaffold for tests, then remove
#repeatMaskCoords=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/repeatMaskingCoordinates/masking_coordinates
outdir=$wd/vcf_filtering
mkdir -p $outdir


##### run through 
# trim alternates and set maxNOCALL to 20% # good
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-L ${mychrom} \
-V ${indir}/${infile} \
-trimAlternates \
--maxNOCALLnumber $maxNoCall \
-o ${outdir}/1_TrimAltRemoveNoCall_${infile}

# testing notes: vcf file requires an index; 

## biallelic snps:
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-L ${mychrom} \
-V ${outdir}/1_TrimAltRemoveNoCall_${infile} \
--restrictAllelesTo BIALLELIC \
--selectTypeToInclude SNP \
-o ${outdir}/'2snp_Filter_TrimAltRemoveNoCall_'${infile}

# testing notes: 
