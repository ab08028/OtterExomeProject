#! /bin/bash
#$ -cwd
#$ -l h_rt=30:00:00,h_data=32G,highp
#$ -N sammpileupBCFcall
#$ -o /u/scratch/a/ab08028/otters/reports
#$ -e /u/scratch/a/ab08028/otters/reports
#$ -m abe
#$ -M ab08028
#$ -t 1-323
# doesn't take very long with samtools! :)
source /u/local/Modules/default/init/modules.sh
module load samtools/1.3.1
module load bcftools/1.3.1
#GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar
# see if it speeds things up to have reference on SCRATCH
# REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/sea_otter_genome/indexGenome/sea_otter_23May2016_bS9RH.fasta
#REFERENCE=/u/scratch2/a/ab08028/otters/sea_otter_23May2016_bS9RH.fasta
# ferret reference!
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fa 
intervalFiles=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/intervalFiles

IN_BAM_DIR=/u/scratch2/a/ab08028/otters/bams/ReadsFiltered
# note needed to separately install tabix:
# and make sure bgzip is in your path. 
tabix=/u/home/a/ab08028/klohmueldata/annabel_data/bin/tabix-0.2.6/tabix

### Alternate step 7: samtools instead of GATK; skipping BQSR
# GATK won't work with de novo genome
# http://www.htslib.org/workflow/#mapping_to_variant
# call variants:
mkdir -p /u/scratch2/a/ab08028/otters/vcfs/Samtools
VCF_DIR=/u/scratch2/a/ab08028/otters/vcfs/Samtools
# going to keep in indel realignment bc samtools doesn't do it:
## also want to keep in monomorphic sites (all callable sites) -- useful down the line. 
samtools mpileup -l $intervalFiles/interval_${SGE_TASK_ID}.bed -Q 20 -E -u -g -f $REFERENCE $IN_BAM_DIR/01_Elut_CA_Gidget_Aligned_MarkDup_IndelRealigned_Filtered.bam | bcftools call -m -O z -o $VCF_DIR/01_Elut_CA_Gidget.raw_variants.${SGE_TASK_ID}.samtools.vcf.gz
# mpileup:
# Q: AB addition to defaults:  -Q 20 , --min-BQ INT
   # Minimum base quality for a base to be considered [13]  
# u: output uncompressed file, good for piping
# g: Compute genotype likelihoods and output them in the binary call format (BCF). 
# f: reference fasta 
# don't use: v: output variant sites only ## take this out. want all sites. 
# E: AB maybe yes (supposed to be better): calculation with -E. IÕve heard that the latter option will be turned on by default in the next version of SAMtools
echo " finished with variant calls; now indexing vcf" 
#
# bcftools call: 
# m: alternative modelfor multiallelic and rare-variant calling designed to overcome known limitations in -c calling model 
# O z: output format z (compressed vcf)
# -o: output file name 
#index vcf
vcf=$VCF_DIR/01_Elut_CA_Gidget.raw_variants.${SGE_TASK_ID}.samtools.vcf.gz
$tabix -p vcf $vcf
# make plots
echo "now making plots" 
bcftools stats -F $REFERENCE $vcf > $vcf.stats
# WITH MORE INDIVIDUALS: -s/-S must given to include also sample columns. with list of samples (not including since only 1 sample)
mkdir $VCF_DIR/plots
plot-vcfstats -p $VCF_DIR/plots/ $vcf.stats
# filter: not doing this yet. 
# bcftools filter -O z -o <study_filtered..vcf.gz> -s LOWQUAL -i'%QUAL>10' <study.vcf.gz>

sleep 5m