# Sandbox : try to get genotype likelihoods with ANGSD 
# Also try with random sampling
# (in the background need to start rerunning paleomix)
captures=/u/flashscratch/a/ab08028/captures/
bamList= # file with 3 bam files 
regionFile= # regions to call genotypes; can get this from GATK define intervals, or just use my target list for elut or my regions bed for mfur
# must be in this format : chr1:1-10000 ; what if there is no chr? Try this on one interval to see if it works.
# Get an interval from the targets file 
testwd=/u/home/a/ab08028/klohmueldata/annabel_data/sandbox-angsd
testRegion="ScbS9RH_100661:10009-11075"
testBamList=$testwd/test.Bam.list
elutRef=/u/home/a/ab08028/klohmueldata/annabel_data/sea_otter_genome/dedup_99_indexed_USETHIS/sea_otter_23May2016_bS9RH.deduped.99.fasta
mfurRef=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta





# covered intervals dir -- merge across all individuals? 

# test Bam list contains pathto : /u/flashscratch/a/ab08028/captures/paleomix/testMapping/original.filter.A30_Elut_CA_SM_35_SN1_CAP/A30_Elut_CA_SM_35_SN1_CAP.sea_otter_23May2016_bS9RH.deduped.99.bam
############### angsd experiment ###########

######################### sea otter ref genome ###########################
# Try angsd:
source activate angsd-conda-env
# ./angsd -GL 2 -trim 4 -nThreads 10 -bam $bamList -rf $regionFile

# testing:
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
-doHaploCall # can either be through random sampling or consensus sampling 

# Could use Arun's ANGSD wrapper
# angsd-wrapper Genotypes ./Genotype_Config

## Unanswered questions
# - what is order of output genotypes? 
# - what are conseqeuences of doMajorMinor 1 vs 4
# - should I use GATK or Samtools style GLs?
# - should i trim 4 or 7bp from reads? Should I do it to modern DNA too for consistency?
# - should I do random haploid sampling? consensus or random?
# - how do I supply ancestral state from ferret? 

# Answered questions:
# - using the Scb scaffold designator works for region format: ScbS9RH_100661:10009-11075 is valid!
# - 

##################### get folded SFS #################
# from angsd docs If you don't have the ancestral state, you can instead estimate the folded SFS. This is done by supplying the -anc with the reference genome and also supply -fold 1. 

#first generate .saf file
angsd -bam $testBamList -trim 4 -nThreads 4 -doSaf 1 -out smallFolded -anc $elutRef -GL 2 -fold 1 -r $testRegion -minQ 20 -minMapQ 30 -skipTriallelic 1 
#now try the EM optimization with 4 threads
realSFS smallFolded.saf.idx -maxIter 100 -P 4 -r $testRegion >smallFolded.sfs
#in R
sfs<-scan("smallFolded.sfs")
barplot(sfs[-1])
######################### ferret ref genome ###########################
