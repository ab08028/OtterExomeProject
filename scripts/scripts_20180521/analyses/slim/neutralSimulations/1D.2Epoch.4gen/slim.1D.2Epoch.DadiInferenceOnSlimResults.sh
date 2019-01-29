#! /bin/bash
#$ -cwd
#$ -l h_rt=02:00:00,h_data=2G
#$ -o /u/flashscratch/a/ab08028/captures/reports/slim
#$ -e /u/flashscratch/a/ab08028/captures/reports/slim
#$ -m abe
#$ -M ab08028
#$ -t 1-100

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
slimModel=1D.2Epoch # model SLIM was simulated under 
gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scripts=$gitdir/scripts/scripts_20180521/analyses
slimscriptdir=$scripts/slim/neutralSimulations/${slimModel}
dadiscriptdir=$scripts/dadi_inference/

rundate= # date slim was run
pop=generic

mu=8.64411385098638e-09
L=6000000 # 6Mb

captures=/u/flashscratch/a/ab08028/captures/
wd=$captures/analyses/slim/neutralSimulations/${slimModel}/${rundate}/replicate_${SGE_TASK_ID}

scripts='1D.1Epoch.dadi.py 1D.2Epoch.dadi.py'


for script in $scripts
do
dadimodel=${script%.dadi.py}
echo "starting inference for $dadimodel"
outdir=$wd/dadiInfBasedOnSlim/dadiInfModel_$dadimodel/
mkdir -p $outdir
# carry out inference with 50 replicates that start with different p0 perturbed params:
for i in {1..50}
do
echo "carrying out inference $i for model $dadimodel for pop $pop" 
# [0-9] indicates that it's a number, but not specific about proj value
python $dadiscriptdir/$script \
--runNum $i \
--pop generic \
--mu $mu \
--L $L \
--sfs $wd/SFS/generic.${slimModel}.slim.output.unfolded.sfs.dadi.format.${rundate}.txt \
--outdir $outdir
done


echo "concatenating results"
grep rundate -m1 $outdir/${pop}.dadi.inference.${model}.runNum.1.*.output > $outdir/${pop}.dadi.inference.${model}.all.output.concatted.txt
for i in {1..50}
do
grep rundate -A1 $outdir/${pop}.dadi.inference.${model}.runNum.${i}.*.output | tail -n1 >> $outdir/${pop}.dadi.inference.${model}.all.output.concatted.txt
done

done

############################ deactivate virtualenv ###############
deactivate # deactivate virtualenv
