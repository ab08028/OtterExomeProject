gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptdir=$gitdir/analyses/slim/

# 100 replicates
for i in {1..100}
do
qsub slim.2Epoch.array.sh $i
done

# this will submit 100 jobs, each of which is an array of 60 portions.