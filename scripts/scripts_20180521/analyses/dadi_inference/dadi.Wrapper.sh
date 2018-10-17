#### wrapper:

scriptdir=/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/analyses/dadi_inference
script=1D.Bottleneck.dadi.py
model="Bottleneck_1D"
mu=8.64411385098638e-09
genotypeDate=20180806
sfsdir=/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/$genotypeDate/neutralSFS/
sfssuffix=all_9_rmAllHet_rmRelativesAdmixed_passingAllFilters_allCalled.filtered.sfs.20181009.dadi.format.out
outdir=/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/$genotypeDate/$model
mkdir -p $outdir
### want to make a slightly fancier outdir that is the model / date or something like that eventually. 
for pop in CA
do
L=4193488 # this is just for CA get L somehow for each population (have to figure this out)
for i in {1..50}
do
python $scriptdir/$script --runNum $i --pop $pop --mu $mu --L $L --sfs ${sfsdir}/${pop}_${sfssuffix} --outdir $outdir
done
done

# set up header
grep rundate -m1 dadi.inference.bottleneck.runNum.1.20181016.output > all.output.concatted.txt
for i in {1..50}
do
grep rundate -A1 dadi.inference.bottleneck.runNum.${i}.20181016.output | tail -n1 >> all.output.concatted.txt
done

# once all are done, you can cat them all together, skipping header for all but first files...?



# test:
