gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
models='1D.2Epoch.1.5Mb.cds.20KAncSize'
populations="genericPop" # do COM separately below
# loop through models, populations and 25 replicates
scriptdir=$gitdir/scripts/scripts_20180521/analyses/slim/cdsSimulations/

todaysdate=`date +%Y%m%d`
for pop in $populations
do
for model in $models
do
# make the slim script from the maker script:
wd=$SCRATCH/captures/analyses/slim/cdsSimulations/$pop/$model/$todaysdate/
mkdir -p $wd
# make a dir to put logs in:
logdir=$wd/logs
mkdir -p $logdir
# make the slim script:
for h in 0 0.5
do
sh $scriptdir/$pop/$model/make_slim_elut.${model}.${pop}.sh $h
done

for i in {1..25}
do
# qsub -N name -o outdir -e errordir $script $pop $model $rep $rundate
# use long-run script:
qsub -N slimRep${i}.${pop} -o $logdir -e $logdir $scriptdir/array_slim_elut.generic.longRun.sh $pop $model $i $todaysdate
done
done
done
