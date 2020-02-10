#! /bin/bash
#$ -cwd
#$ -l h_rt=24:00:00,h_data=8G
#$ -N dadi_inference
#$ -o /u/flashscratch/a/ab08028/captures/reports/dadi
#$ -e /u/flashscratch/a/ab08028/captures/reports/dadi
#$ -m abe
#$ -M ab08028


######################### set up virtual environment #####################
# if using on Hoffman, have to use virtual environment:
# you can set up the virtual env by doing:
# module load python/2.7.13_shared
# virtualenv $HOME/env_python2.7.13 # once
# then activate it every future time with  source $HOME/env_python2.7.13/bin/activate
# make sure you pip install numpy scipy matplotlib 
# and then you can deactivate

source /u/local/Modules/default/init/modules.sh
module load python/2.7.13_shared
source /u/home/a/ab08028/env_python2.7.13/bin/activate # activate virtual environment 

SCRATCH=/u/flashscratch/a/ab08028/
gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/ # hoffman
scriptdir=$gitdir/scripts/scripts_20180521/analyses/dadi_inference
#script=$1 # generic -- give it any of your dadi scripts -- make a version of this for Hoffman too.
#model=${script%.dadi.py}
mu=8.64411385098638e-09
genotypeDate=20181119 # newer gts
sfsDate=20181221 # projection with 0.75 het filter and these projection values
hetFilter=0.75

todaysdate=20200210

captures=$SCRATCH/captures/
sfsdir=$captures/analyses/SFS/$genotypeDate/easySFS/neutral/projection-${sfsDate}-hetFilter-${hetFilter}/dadi-plusMonomorphic/
dadidir=$captures/analyses/dadi_inference/
sfssuffix=plusMonomorphic.sfs
### Make sure this is the correct file #####
#totalNeut=$captures/vcf_filtering/${genotypeDate}_filtered/bedCoords/neutralCallableSites_perPop/summary.neutralCallableSites.perPop.txt # file with total neutral sites counts for each population 
### want to make a slightly fancier outdir that is the model / date or something like that eventually. 
# run multiple models for multiple popuations?
#scripts='1D.1Bottleneck.dadi.py 1D.2Bottleneck.dadi.py 1D.2Epoch.dadi.py' # list of models you want to run
#scripts='1D.2Epoch.dadi.py 1D.1Bottleneck.TB20gen.dadi.py'
# want to run with new starting parameters
#scripts='1D.2Epoch.newStartValues.dadi.py 1D.2Epoch.LargeFoldChange.dadi.py'
#scripts='1D.1Epoch.dadi.py' # just this one for now
#for pop in CA AK AL COM KUR
#scripts='1D.1Epoch.dadi.py' # just this one for now
#scripts='2D.Bottleneck.Migration.dadi.py 2D.ConstantSize.Migration.dadi.py 2D.ConstantSize.noMigration.dadi.py'
scripts='2D.ConstantSize.noMigration.dadi.py'
for popPair in CA-AK CA-KUR AK-KUR AK-AL AL-KUR AK-COM AL-COM CA-AL CA-COM CA-KUR # 
do
echo $popPair
# get total sites from total sites file that was written out as part of my easySFS scripts
L=`grep -v totalSites $sfsdir/${popPair}.totalSiteCount.L.withMonomorphic.txt | awk '{print $3}'` # get combined total sites

for script in $scripts
do

model=${script%.dadi.py}
echo "starting inference for $popPair for model $model"
outdir=$dadidir/$genotypeDate/$popPair/inference_$todaysdate/$model/
mkdir -p $outdir
# carry out inference with 50 replicates that start with different p0 perturbed params:
for i in {1..50}
do
echo "carrying out inference $i for model $model for pop $popPair" 
# [0-9] indicates that it's a number, but not specific about proj value
python $scriptdir/$script --runNum $i --pop $popPair --mu $mu --L $L --sfs ${sfsdir}/${popPair}.${sfssuffix} --outdir $outdir
done


echo "concatenating results"
grep rundate -m1 $outdir/${popPair}.dadi.inference.${model}.runNum.1.output > $outdir/${popPair}.dadi.inference.${model}.all.output.concatted.txt
for i in {1..50}
do
grep rundate -A1 $outdir/${popPair}.dadi.inference.${model}.runNum.${i}.output | tail -n1 >> $outdir/${popPair}.dadi.inference.${model}.all.output.concatted.txt
done

done
done
############################ deactivate virtualenv ###############
deactivate # deactivate virtualenv
