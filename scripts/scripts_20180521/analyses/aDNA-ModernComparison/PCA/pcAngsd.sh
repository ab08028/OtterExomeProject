#! /bin/bash
#$ -cwd
#$ -l h_rt=100:00:00,h_data=2G,highp
#$ -m abe
#$ -M ab08028
#$ -pe shared 10
#$ -N PCAngsd
#$ -e /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -o /u/flashscratch/a/ab08028/captures/reports/angsd

##### pcangsd -- run PCAngsd to do PCA based on genotype likelihoods
# set a MAF of 0.05
# maybe do some sort of SNP likelihood filtering? 1e-6?

source /u/local/Modules/default/init/modules.sh
module load anaconda # load anaconda
source activate angsd-conda-env # activate conda env
module load python
pcangsddir=/u/home/a/ab08028/klohmueldata/annabel_data/bin/pcangsd
wd=$SCRATCH/captures/aDNA-ModernComparison
angsdDate=20190502 # date GLs were called
#GLdir=$wd/angsd-GLs/$angsdDate # eventually
# for now
GLdir=/u/flashscratch/a/ab08028/captures/aDNA-ModernComparison/angsd-GLs #temporary!
PCAdir=$wd/pca/covarianceMatrices
mkdir -p $PCAdir
for ref in Elut Mfur
do
python $pcangsddir/pcangsd.py \
-beagle $GLdir/angsdOut.mappedToElut.beagle.gz \
-o $PCAdir/pcAngsd.$ref \
-minMaf 0.05 \
-threads 10
# default minMaf is 0.05

done


# then have to move to python or R to deal with this
# using a modification of:
# https://github.com/mfumagalli/ngsPopGen/blob/master/scripts/plotPCA.R