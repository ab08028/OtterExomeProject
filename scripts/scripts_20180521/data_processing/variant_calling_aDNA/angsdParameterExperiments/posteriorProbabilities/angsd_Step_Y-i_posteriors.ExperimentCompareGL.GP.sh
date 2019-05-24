#! /bin/bash
#$ -cwd
#$ -l h_rt=200:00:00,h_data=2G,highp
#$ -m abe
#$ -M ab08028
#$ -pe shared 16
#$ -N angsdOrlandoFullCov
#$ -e /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -o /u/flashscratch/a/ab08028/captures/reports/angsd

# idea:
# pick a chromosome of ferret: GL896898.1
# and output beagle GLs and GPs 
# based on high cov and ancient only
# then separately with low cov and ancient only
# going to plot GL vs GP for each Genotype separately (Arun suggestion)
#### ANGSD v 0.923 ####
source /u/local/Modules/default/init/modules.sh
module load anaconda # load anaconda
source activate angsd-conda-env # activate conda env
gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptDir=$gitDir/scripts/scripts_20180521/data_processing/variant_calling_aDNA
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/aDNA-ModernComparison
bamdir=$wd/bams/
GLdir=$wd/angsd-GLs
mkdir -p $GLdir
#todaysdate=`date +%Y%m%d`
################### high coverage #########################
todaysdate="20190523-experiment-GLvsGP-highcov"
mkdir -p $GLdir/$todaysdate/posteriorProbabilities
outdir=$GLdir/$todaysdate/posteriorProbabilities

scaff="GL896898.1:1-3000000" # only first 3M sites -- keeping it small for now 


mfurBamList=$scriptDir/bamLists/angsd.bamList.mappedtoMfurfullpaths.HighCovPlusADNAOnly.txt  # list of bam files mapped to ferret, including downsampled AND non-downsampled

# references:
mfurRef=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta

echo "THIS USES HIGH COVERAGE MODERN + ANCIENT ONLY, only $scaff" > $GLdir/$todaysdate/HIGHCOVERAGEONLY.$scaff.txt
# trying output in beagle format  doGlf 2
####### Mfur mapped bams ############
spp="mfur"
ref=$mfurRef
bamList=$mfurBamList
# doGlf 2 to output beagle format GLs
angsd -nThreads 16 \
-ref $ref \
-bam $bamList \
-GL 2 \
-doMajorMinor 1 -doMaf 1 \
-beagleProb 1 -doPost 1 \
-remove_bads 1 -uniqueOnly 1 \
-C 50 -baq 1 -trim 4 -minQ 20 -minMapQ 25 -skipTriallelic 1 \
-out $outdir/angsdOut.mappedTo${spp}.OrlandoSettings \
-doGlf 2 \
-r $scaff \
-doCounts 1 -dumpCounts 2
# 750,800 sites
# adding detail that it outputs post-filtering counts per site per individual (doCounts 1 dumpCounts 2)

################### low coverage #########################
todaysdate="20190523-experiment-GLvsGP-lowcov"
mkdir -p $GLdir/$todaysdate/posteriorProbabilities
outdir=$GLdir/$todaysdate/posteriorProbabilities

scaff="GL896898.1:1-3000000"


mfurBamList=$scriptDir/bamLists/angsd.bamList.LowCoverageOnly.mappedtoMfurfullpaths.txt # low coverage + aDNA

# references:
mfurRef=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta

echo "THIS USES DOWNSAMPLED LOW COVERAGE MODERN + ANCIENT ONLY, only $scaff" > $GLdir/$todaysdate/LOWCOVERAGEONLY.$scaff.txt
# trying output in beagle format  doGlf 2
####### Mfur mapped bams ############
spp="mfur"
ref=$mfurRef
bamList=$mfurBamList
# doGlf 2 to output beagle format GLs
angsd -nThreads 16 \
-ref $ref \
-bam $bamList \
-GL 2 \
-doMajorMinor 1 -doMaf 1 \
-beagleProb 1 -doPost 1 \
-remove_bads 1 -uniqueOnly 1 \
-C 50 -baq 1 -trim 4 -minQ 20 -minMapQ 25 -skipTriallelic 1 \
-out $outdir/angsdOut.mappedTo${spp}.OrlandoSettings \
-doGlf 2 \
-r $scaff \
-doCounts 1 -dumpCounts 2

# 99,711 sites

######## then going to take the GLs and GPs and plot in RR like a QQ plot
### will be interesting
# do for 0/0 0/1 and 1/1
# for all individuals
# yes.
# maybe also bring in the mafs
# and plot maf vs GP or something?
