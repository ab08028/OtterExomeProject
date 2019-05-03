#! /bin/bash
#$ -cwd
#$ -l h_rt=100:00:00,h_data=2G,highp
#$ -m abe
#$ -M ab08028
#$ -pe shared 16
#$ -N angsdGLs
#$ -e /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -o /u/flashscratch/a/ab08028/captures/reports/angsd

#### ANGSD ####
source /u/local/Modules/default/init/modules.sh
module load anaconda # load anaconda
source activate angsd-conda-env # activate conda env
gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptDir=$gitDir/scripts/scripts_20180521/data_processing/variant_calling_aDNA/
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/aDNA-ModernComparison
bamdir=$wd/bams/
GLdir=$wd/angsd-GLs
mkdir -p $GLdir
todaysdate=`date +%Y%m%d`
mkdir -p $GLdir/$todaysdate
# this is temporary -- just calling in one region to make sure angsd works
# then maybe want to call genome-wide whereever we can?
# or restrict to called regions 
testRegion="ScbS9RH_100661:10009-11075"

# gather bams from paleomix using script gatherBamsForDownsampling.sh
# and make lists of the relevant bams: 

elutBamList=$scriptDir/bamLists/angsd.bamList.mappedtoElutfullpaths.txt # list of bam files mapped to sea otter, including downsampled AND non-downsampled
mfurBamList=$scriptDir/bamLists/angsd.bamList.mappedtoMfurfullpaths.txt  # list of bam files mapped to ferret, including downsampled AND non-downsampled

# references:
elutRef=/u/home/a/ab08028/klohmueldata/annabel_data/sea_otter_genome/dedup_99_indexed_USETHIS/sea_otter_23May2016_bS9RH.deduped.99.fasta
mfurRef=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta

# trying output in beagle format  doGlf 2
####### Mfur mapped bams ############
angsd \
-GL 2 \
-trim 4 \
-nThreads 16 \
-bam $mfurBamList \
-minQ 20 -minMapQ 30 \
-skipTriallelic 1 \
-doMajorMinor 4 -ref $mfurRef \
-doGlf 2 \
-uniqueOnly 1 \
-doMaf 2 \
-out $GLdir/$todaysdate/angsdOut.mappedToMfur.allSites \
-remove_bads 1 \
-C 50

# 20190502 -- was run without doDepth or doCount
# 201090503 -- going to add more things: doDepth/doCount to get depth per sample
# Adding more filtering :
# remove_bads and -C50

########### Elut mapped bams #####################
angsd \
-GL 2 \
-trim 4 \
-nThreads 16 \
-bam $elutBamList \
-minQ 20 -minMapQ 30 \
-skipTriallelic 1 \
-doMajorMinor 4 -ref $elutRef \
-doGlf 2 \
-uniqueOnly 1 \
-doMaf 2 \
-out $GLdir/$todaysdate/angsdOut.mappedToElut.allSites \
-remove_bads 1 \
-C 50
# not sure: -only_proper_pairs if I should use or not... 

# 20190502 -- was run without doDepth or doCount
# 201090503 -- going to add more things: doDepth/doCount to get depth per sample
# Adding more filtering :
# remove_bads and -C50


source deactivate

sleep 10m

############ info on flags: ##########
# doMaf 2 -- fixed major and unknown minor  "Known major, Unknown minor. Here the major allele is assumed to be known (inferred or given by user) however the minor allele is not determined. Instead we sum over the 3 possible minor alleles weighted by their probabilities. T"
#uniqueOnly -- removes reads with more than 1 best hit
# doGlf is how to do output. 4 is gzipped text
# -doMajorMinor 4 -ref $elutRef sets the major and minor alleles based on the reference (like GATK); makes less sense for elut-mapped bc california isn't ancestral. makes more sense for ferret
# Major Minor is weird You can force the major allele according to the reference states if you have defined those -ref. The minor allele will be inferred based on the genotype likelihood (see do major minor 1). This is the approach used by both GATK and Samtools
# what are consequences
#-doMajorMinor 4
#-ref [fasta.fa]
# more from Daly 
# nThreads : number of threads
# trying with GL2 =  GATK style GLs for now, but also want to try with samtools (GL 1)
# trim =4 ; trim 4bp on either end of reads (to get rid of damage)
# output is chr, position, then the LogLikelihoods for each individual for each possible genotype AA AC AG AT CC CG CT GG GT TT
### CAN also do haploid calls by sampling a single read or taking consensus (what is best approach?)
# -doHaploCall # can either be through random sampling or consensus sampling 

## Unanswered questions
# - what is order of output genotypes? 
# - what are conseqeuences of doMajorMinor 1 vs 4
# - should I use GATK or Samtools style GLs?
# - should i trim 4 or 7bp from reads? Should I do it to modern DNA too for consistency?
# - should I do random haploid sampling? consensus or random? <-- I want to do this too maybe.
# - how do I supply ancestral state from ferret? 
