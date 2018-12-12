# Gather up FSC results

# want a table for each model/population with all results
models='1D.1Bottleneck'
pops="CA AK AL COM KUR"
genotypeDate=20180806 # date genotypes were called
sfsDate=20181019 # date sfses were made
rundate=20181031 # date you are running inference (to distinguish later analyses) <-- can set this with date command or manually

gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptDir=$gitDir/scripts/scripts_20180521/analyses/
wd=/u/flashscratch/a/ab08028/captures/analyses/ # overall working dir
infDir=$wd/fastsimcoal_inference # specify where inference is happening

mkdir $infDir/resultsSummaries

for pop in $pops
do
############ get sample size for the population ############
sumdir=$infDir/resultsSummaries/$pop
mkdir -p $sumdir

for model in $models
do
header=${model}

outfile=$sumdir/${pop}.fsc.inference.${model}.${rundate}.all.output.concatted.txt
# get header:
header=`head -n1 $infDir/$pop/inference_${rundate}/$model/run_1/$model/*bestlhoods`
echo -e "runNum\t$header" > $outfile
########### copy generic files into directory and update #########
for i in {1..50}
do
outdir=$infDir/$pop/inference_${rundate}/$model/run_${i}/$model # specify where output will be pulled from
results=`grep -v [A-Z] $outdir/*bestlhoods`
echo -e "$i\t$results" >> $outfile


# try

done
done
done
