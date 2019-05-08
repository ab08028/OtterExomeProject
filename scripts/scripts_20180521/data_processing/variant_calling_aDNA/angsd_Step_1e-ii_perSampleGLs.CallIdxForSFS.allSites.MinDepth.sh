#! /bin/bash
#$ -cwd
#$ -l h_rt=10:00:00,h_data=2G,highp
#$ -m abe
#$ -M ab08028
#$ -pe shared 8
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
mkdir -p $SFSdir/$todaysdate/perIndividual-minDepth
mkdir -p $GLdir/$todaysdate/perIndividual-minDepth

######## info from command line : # ##
bam=$1
label=$2 # ID.ref.downsample.rep#
refPrefix=$3
reference=$4
minDepth=$5
maxDepth=$6
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
-out $GLdir/$todaysdate/perIndividual-minDepth/${label}.mappedTo${refPrefix}.allSites.minDepth.$minDepth.maxDepth.$maxDepth \
-remove_bads 1 \
-C 50 \
-doCount \
-SetMinDepth $minDepth \
-SetMaxDepth $maxDepth


# removed : -SNP_pval 1e-6 \

# get it without transitions: 
angsd \
-dosaf 1 \
-fold 1 \
-noTrans 1 \
-anc $reference \
-ref $reference \
-fai ${reference}.fai \
-glf $GLdir/$todaysdate/perIndividual-minDepth/${label}.mappedTo${refPrefix}.allSites.minDepth.$minDepth.maxDepth.$maxDepth.glf.gz \
-out $SFSdir/$todaysdate/perIndividual-minDepth/${label}.mappedTo${refPrefix}.allSites.minDepth.$minDepth.maxDepth.$maxDepth.TransversionsOnly \
-nInd 1

# get it with transitions+transversions (better for pi estimate?): 
angsd \
-dosaf 1 \
-fold 1 \
-anc $reference \
-ref $reference \
-fai ${reference}.fai \
-glf $GLdir/$todaysdate/perIndividual-minDepth/${label}.mappedTo${refPrefix}.allSites.minDepth.$minDepth.maxDepth.$maxDepth.glf.gz \
-out $SFSdir/$todaysdate/perIndividual-minDepth/${label}.mappedTo${refPrefix}.allSites.minDepth.$minDepthmaxDepth.$maxDepth \
-nInd 1

realSFS $SFSdir/$todaysdate/perIndividual-minDepth/${label}.mappedTo${refPrefix}.allSites.minDepth.$minDepth.maxDepth.$maxDepth.TransversionsOnly.saf.idx > $SFSdir/$todaysdate/perIndividual-minDepth/${label}.mappedTo${refPrefix}.allSites.minDepth.$minDepth.maxDepth.$maxDepth.TransversionsOnly.saf.SFS.txt
realSFS $SFSdir/$todaysdate/perIndividual-minDepth/${label}.mappedTo${refPrefix}.allSites.minDepth.$minDepth.maxDepth.$maxDepth.saf.idx > $SFSdir/$todaysdate/perIndividual-minDepth/${label}.mappedTo${refPrefix}.allSites.minDepth.$minDepth.maxDepth.$maxDepth.saf.SFS.txt

source deactivate

sleep 10m