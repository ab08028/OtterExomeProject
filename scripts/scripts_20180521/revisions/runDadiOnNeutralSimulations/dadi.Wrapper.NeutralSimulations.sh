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
#sfsDate=20190301 # projection with 0.75 het filter and these projection values
#hetFilter=0.75
todaysdate=`date +%Y%m%d`
captures=$SCRATCH/captures/
#sfsdir=$captures/analyses/SFS/$genotypeDate/easySFS/neutral/projection-${sfsDate}-hetFilter-${hetFilter}/dadi-plusMonomorphic/
sfsdir=/u/home/a/ab08028/klohmueldata/annabel_data/captures/analyses/slim/neutralSimulations_containsTarballs/neutralSimsForManuscript_20200511/CA.1D.2Epoch.35Gen.200Inds/20200224/allSFSes
dadidir=$captures/analyses/dadi_inference/
#sfssuffix=plusMonomorphic.sfs
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
scripts='1D.1Epoch.dadi.py 1D.2Epoch.dadi.py'
simulationModel='CA.1D.2Epoch.35Gen.200Inds.postContraction'
for pop in CA # skipping med/ber/com for now
do
# get total sites from total sites file that was written out as part of my easySFS scripts
#L=`grep $pop $sfsdir/$pop-[0-9]*.totalSiteCount.L.withMonomorphic.txt | awk '{print $2}'`
L=6000000 # simulated 6Mb
for script in $scripts
do
for rep in 1 2 3 4 5 6 7 8 9 11
do
model=${script%.dadi.py}
echo "starting inference for $pop for model $model"
outdir=$dadidir/runningDadiOnNeutralSimulations/$simulationModel/inference_$todaysdate/$model/replicate_${rep}
mkdir -p $outdir
# carry out inference with 50 replicates that start with different p0 perturbed params:
for i in {1..50}
do
echo "carrying out inference $i for model $model for pop $pop" 
# [0-9] indicates that it's a number, but not specific about proj value
python $scriptdir/$script --runNum $i --pop $pop --mu $mu --L $L --sfs $sfsdir/$pop.replicate_${rep}.${simulationModel}.slim.output.unfolded.sfs.dadi.format.txt --outdir $outdir
done


echo "concatenating results"
grep rundate -m1 $outdir/${pop}.dadi.inference.${model}.runNum.1.output > $outdir/${pop}.dadi.inference.${model}.all.output.concatted.txt
for i in {1..50}
do
grep rundate -A1 $outdir/${pop}.dadi.inference.${model}.runNum.${i}.output | tail -n1 >> $outdir/${pop}.dadi.inference.${model}.all.output.concatted.txt
done

done
done
done
############################ deactivate virtualenv ###############
deactivate # deactivate virtualenv
