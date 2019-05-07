#! /bin/bash
#$ -cwd
#$ -l h_rt=4:00:00,h_data=2G
#$ -m abe
#$ -M ab08028
#$ -pe shared 4
#$ -N angsdSFSPerInd
#$ -e /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -o /u/flashscratch/a/ab08028/captures/reports/angsd

#### ANGSD v 0.923 ####
source /u/local/Modules/default/init/modules.sh
module load anaconda # load anaconda
source activate angsd-conda-env # activate conda env
gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptDir=$gitDir/scripts/scripts_20180521/data_processing/variant_calling_aDNA/
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/aDNA-ModernComparison
bamdir=$wd/bams/
GLdir=$wd/angsd-GLs
SFSdir=$wd/angsd-SFS
mkdir -p $SFSdir
mkdir -p $GLdir
todaysdate=`date +%Y%m%d`
snpCutoff=1e-6
mkdir -p $GLdir/$todaysdate
mkdir -p $SFSdir/$todaysdate/perIndividual
mkdir -p $GLdir/$todaysdate/perIndividual

######## info from command line : # ##

bam=$1
sampleID=$2
refPrefix=$3
reference=$4
### basename pulls the file name out of the bath so you can use that 
angsd \
-GL 2 \
-trim 4 \
-nThreads 8 \
-bam $bam \
-minQ 20 -minMapQ 30 \
-skipTriallelic 1 \
-doMajorMinor 4 -ref $reference \
-doGlf 1 \
-uniqueOnly 1 \
-doMaf 2 \
-out $GLdir/$todaysdate/perIndividual/${bam%.bam}.mappedTo${refPrefix}.allSites \
-remove_bads 1 \
-C 50 
# removed : -SNP_pval 1e-6 \

# get it without transitions: 
angsd \
-dosaf 1 \
-fold 1 \
-noTrans 1 \
-anc $reference \
-ref $reference \
-fai ${reference}.fai \
-glf $GLdir/$todaysdate/perPopulation/${bam%.bam}.mappedTo${refPrefix}.allSites.glf.gz \
-out $SFSdir/$todaysdate/${bam%.bam}.mappedTo${refPrefix}.allSites.TransversionsOnly \
-nInd 1

# get it with transitions+transversions (better for pi estimate?): 
angsd \
-dosaf 1 \
-fold 1 \
-anc $reference \
-ref $reference \
-fai ${reference}.fai \
-glf $GLdir/$todaysdate/perPopulation/${bam%.bam}.mappedTo${refPrefix}.allSites.glf.gz \
-out $SFSdir/$todaysdate/${bam%.bam}.mappedTo${refPrefix}.allSites \
-nInd 1

realSFS $SFSdir/$todaysdate/${bam%.bam}.mappedTo${refPrefix}.allSites.TransversionsOnly.saf.idx > $SFSdir/$todaysdate/perIndividual/${bam%.bam}.mappedTo${refPrefix}.allSites.TransversionsOnly.saf.SFS.txt
realSFS $SFSdir/$todaysdate/${bam%.bam}.mappedTo${refPrefix}.allSites.saf.idx > $SFSdir/$todaysdate/perIndividual/${bam%.bam}.mappedTo${refPrefix}.allSites.saf.SFS.txt

done

source deactivate

sleep 10m
