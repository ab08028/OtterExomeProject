#! /bin/bash
#$ -cwd
#$ -l h_rt=03:00:00,h_data=8G
#$ -o /u/flashscratch/a/ab08028/captures/reports/slim
#$ -e /u/flashscratch/a/ab08028/captures/reports/slim
#$ -m abe
#$ -M ab08028
#$ -t 1-11
#$ -N dadiInfOnSlim

# 100 replicates (or however many you did with slim)
## order
# slim: 100 replicates of 6Mb (made up 100 chunks, that are concatenated into one output file) --> SFS (one per replicate)
# --> dadi inference based on SFS (1Epoch and 2Epoch) --> 50 dadi replicates.
# so it's 50 dadi replicates per slim replicate , so 50 *100 = 5000 total dadi runs per model (100,000 across two models)


######### do dadi inference based on replicates ##########
source /u/local/Modules/default/init/modules.sh
module load python/2.7.13_shared
source /u/home/a/ab08028/env_python2.7.13/bin/activate # activate virtual environment 


# do 2epoch and 1 epoch inference on the 2epoch model
slimModel=1D.1Epoch # model SLIM was simulated under 
gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scripts=$gitdir/scripts/scripts_20180521/analyses
slimscriptdir=$scripts/slim/neutralSimulations/${slimModel}
dadiscriptdir=$scripts/dadi_inference/

rundate=20190125 # date slim was run
pop=generic

mu=8.64411385098638e-09
L=6000000 # 6Mb

captures=/u/flashscratch/a/ab08028/captures/
wd=$captures/analyses/slim/neutralSimulations/${slimModel}/${rundate}/
#repdir=$wd/replicate_${SGE_TASK_ID}
sfsdir=$wd/allSFSes
mkdir -p $wd/dadiInfBasedOnSlim/allDadiResultsConcatted
scripts='1D.1Epoch.dadi.py 1D.2Epoch.dadi.py'


for script in $scripts
do
dadimodel=${script%.dadi.py}
echo "starting inference for $dadimodel"
outdir=$wd/dadiInfBasedOnSlim/dadiInfModel_$dadimodel/replicate_${SGE_TASK_ID}
concatdir=$wd/dadiInfBasedOnSlim/allDadiResultsConcatted/dadiInfModel_$dadimodel/
mkdir -p $outdir
mkdir -p $concatdir
# carry out inference with 50 replicates that start with different p0 perturbed params:
for i in {1..50}
do
echo "carrying out inference $i for model $dadimodel for pop $pop" 
# [0-9] indicates that it's a number, but not specific about proj value
python $dadiscriptdir/$script \
--runNum $i \
--pop $pop \
--mu $mu \
--L $L \
--sfs $sfsdir/${pop}.rep.${SGE_TASK_ID}.${slimModel}.slim.output.unfolded.sfs.dadi.format.txt \
--outdir $outdir
done
# note the date on the sfs is the date it was made ; not super helpful. can I fix that in the python script?


echo "concatenating results"
grep rundate -m1 $outdir/${pop}.dadi.inference.${dadimodel}.runNum.1.*.output > $outdir/${pop}.rep.${SGE_TASK_ID}.dadi.inf.${dadimodel}.all.output.concatted.txt
for i in {1..50}
do
grep rundate -A1 $outdir/${pop}.dadi.inference.${dadimodel}.runNum.${i}.*.output | tail -n1 >> $outdir/${pop}.rep.${SGE_TASK_ID}.dadi.inf.${dadimodel}.all.output.concatted.txt
done

# copy results to the concat folder for download
/bin/cp -f $outdir/${pop}.rep.${SGE_TASK_ID}.dadi.inf.${dadimodel}.all.output.concatted.txt $concatdir/

done

############################ deactivate virtualenv ###############
deactivate # deactivate virtualenv

sleep 10m
