model=1D.2Epoch.CA
gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptdir=$gitdir/scripts/scripts_20180521/analyses/slim/neutralSimulations/populationSpecificModels/${model}
todaysdate=`date +%Y%m%d` # don't want to use todays date because then different time starting arrays could get messed up
# 100 replicates; starting with 10
for i in {1..11}
do
# name job slimRepX
qsub -N slimRep$i $scriptdir/slim.${model}.array.sh $i $todaysdate
done

# this will submit 100 jobs, each of which is an array of 60 portions.

