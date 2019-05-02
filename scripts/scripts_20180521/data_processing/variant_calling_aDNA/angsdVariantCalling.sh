#### ANGSD ####

source activate angsd-conda-env # activate conda env


elutBamList= # list of bam files mapped to sea otter, including downsampled AND non-downsampled
mfurBamList= # list of bam files mapped to ferret, including downsampled AND non-downsampled
elutRef=/u/home/a/ab08028/klohmueldata/annabel_data/sea_otter_genome/dedup_99_indexed_USETHIS/sea_otter_23May2016_bS9RH.deduped.99.fasta
mfurRef=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta

#### 		####

# make lists # 
########### Elut mapped bams #####################
angsd -GL 2 -trim 4 -nThreads 4 -bam $testBamList -r $testRegion -minQ 20 -minMapQ 30 -skipTriallelic 1 -doMajorMinor 4 -ref $elutRef -doGlf 4
# -doMajorMinor 4 -ref $elutRef doesn't make that much sense for elut; does make sense for ferret; 
# Major Minor is weird You can force the major allele according to the reference states if you have defined those -ref. The minor allele will be inferred based on the genotype likelihood (see do major minor 1). This is the approach used by both GATK and Samtools
# what are consequences
#-doMajorMinor 4
#-ref [fasta.fa]
# omg this worked! what does output mean?
# more from Daly 
# nThreads : number of threads
# trying with GL2 =  GATK style GLs for now, but also want to try with samtools (GL 1)
# trim =4 ; trim 4bp on either end of reads (to get rid of damage)
# output is chr, position, then the LogLikelihoods for each individual for each possible genotype AA AC AG AT CC CG CT GG GT TT
### CAN also do haploid calls by sampling a single read or taking consensus (what is best approach?)
# -doHaploCall # can either be through random sampling or consensus sampling 

# Could use Arun's ANGSD wrapper
# angsd-wrapper Genotypes ./Genotype_Config

## Unanswered questions
# - what is order of output genotypes? 
# - what are conseqeuences of doMajorMinor 1 vs 4
# - should I use GATK or Samtools style GLs?
# - should i trim 4 or 7bp from reads? Should I do it to modern DNA too for consistency?
# - should I do random haploid sampling? consensus or random?
# - how do I supply ancestral state from ferret? 