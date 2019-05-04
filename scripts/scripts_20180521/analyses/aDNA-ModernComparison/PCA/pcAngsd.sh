#! /bin/bash
#$ -cwd
#$ -l h_rt=10:00:00,h_data=2G,highp
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
pcangsddir=/u/home/a/ab08028/klohmueldata/annabel_data/bin/pcangsd
wd=$SCRATCH/captures/aDNA-ModernComparison
angsdDate=20190503 # date GLs were called
GLdir=$wd/angsd-GLs/$angsdDate # eventually
# for now
#GLdir=/u/flashscratch/a/ab08028/captures/aDNA-ModernComparison/angsd-GLs #temporary!
PCAdir=$wd/pca/covarianceMatrices/$angsdDate
mkdir -p $PCAdir
for state in 1e-6.snpsOnly 1e-6.snpsOnly.transversionsOnly 1e-6.snpsOnly.downSampOnly.minInd.5 1e-6.snpsOnly.downSampOnly.minInd.5.transversionsOnly  # allSites don't use allSites for PCA, just use SNPs. use with and without transversions
do
for ref in Elut Mfur

do
## is there something I can do here to get minIndividuals? maybe has to be from angsd itself 
angsd \
-beagle $GLdir/angsdOut.mappedTo${ref}.${state}.beagle.gz \
-minInd 
python $pcangsddir/pcangsd.py \
-beagle $GLdir/angsdOut.mappedTo${ref}.${state}.beagle.gz \
-o $PCAdir/pcAngsd.$ref.$state \
-minMaf 0.05 \
-threads 10

## want to do some sort of missingness filter -- maybe sites that are called in all individuals? or in 13/15 or something, so it can't
# be all ancient vs all modern? 
# default minMaf is 0.05
# this generates a covariance matrix called pcAngsd.ref.state.cov.npy which is a numpy binary file
done
done



source deactivate

# then have to move to python or R to deal with this
# using a modification of:
# https://github.com/mfumagalli/ngsPopGen/blob/master/scripts/plotPCA.R


# then can
# using a modification of:
# https://github.com/mfumagalli/ngsPopGen/blob/master/scripts/plotPCA.R
# to plot the pca
# using script plot.pcAngsd.R
