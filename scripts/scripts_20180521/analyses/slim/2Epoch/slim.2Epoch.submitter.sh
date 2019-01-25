model=2Epoch
gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptdir=$gitdir/scripts/scripts_20180521/analyses/slim/$model
todaysdate=`date +%Y%m%d` # don't want to use todays date because then different time starting arrays could get messed up
# 100 replicates
for i in {1..100}
do
# name job slimRepX
qsub -N slimRep$i $scriptdir/slim.2Epoch.array.sh $i $todaysdate
done

# this will submit 100 jobs, each of which is an array of 60 portions.

