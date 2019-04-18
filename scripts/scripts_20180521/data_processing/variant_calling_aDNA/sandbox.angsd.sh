# Sandbox : try to get genotype likelihoods with ANGSD 
# Also try with random sampling
# (in the background need to start rerunning paleomix)

bamList= # file with 3 bam files 
region= # regions to call genotypes; can get this from GATK define intervals, or just use my target list for elut or my regions bed for mfur
# must be in this format : chr1:1-10000 ; what if there is no chr? Try this on one interval to see if it works.
# Get an interval from the targets file 


############### angsd experiment ###########

######################### sea otter ref genome ###########################
# Try angsd:
source activate angsd-conda-env
./angsd -GL 2 -trim 4 -nThreads 10 -bam $bamList -rf $region
# nThreads : number of threads
# trying with GL2 =  GATK style GLs for now, but also want to try with samtools (GL 1)
# trim =4 ; trim 4bp on either end of reads (to get rid of damage)



# Could use Arun's ANGSD wrapper
angsd-wrapper Genotypes ./Genotype_Config



######################### ferret ref genome ###########################
