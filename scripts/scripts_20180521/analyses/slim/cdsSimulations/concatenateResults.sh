# process output of slim 
gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
models='1D.2Epoch.1.5Mb.cds'
# choose the specific combination of populations/models/rundates you want? this is awkward... what is best way to do it?
modelsRundates='1D.2Epoch.1.5Mb.cds/20190404' # maybe? -- this is kind of awkward, maybe have to deal with diff populations differently?
#populations='AK AL CA COM KUR'
#populations="AK AL CA KUR" # do COM separately below
populations="AK"
# loop through models, populations and 25 replicates
scriptdir=$gitdir/scripts/scripts_20180521/analyses/slim/cdsSimulations/

rundate=`date +%Y%m%d`

for pop in $populations
do
echo $pop
for modelRundate in $modelsRundates
do
# make the slim script from the maker script:
wd=$SCRATCH/captures/analyses/slim/cdsSimulations/ # combine model and rundate
vcfOutDir=$wd/concattedVCFs/$pop/$modelRundate
summaryOutDir=$wd/concattedSummaries/$pop/$modelRundate
mkdir -p $vcfOutDir
mkdir -p $summaryOutDir
for state in Pre Post # pre and post contraction
do
for i in {1..1}
do
echo "rep $i ${state}Contraction"
repdir=$SCRATCH/captures/analyses/slim/cdsSimulations/$pop/$modelRundate/replicate_${i}

outSummary=$summaryOutDir/rep.${i}.slim.output.${state}Contraction.allConcatted.summary.txt

outVCF=$vcfOutDir/rep.${i}.slim.output.${state}Contraction.allConcatted.vcf
#vcf header:
grep "#" $repdir/slim.output.${state}Contraction.1.vcf > $outVCF

# summary header:
grep "replicate" $repdir/slim.output.${state}Contraction.1.summary.txt > $outSummary

# then loop over all chunks
for j in {1..20}
do
# concat summaries and gzip
#echo "concatenating summaries"
#  grep -v "^$" removes extra blank lines
grep -v "replicate" $repdir/slim.output.${state}Contraction.${j}.summary.txt | grep -v "^$" >> $outSummary
# concat vcfs and gzip
#echo "concatenating vcf chunks"
# # want to select everything but first column (setting it to "") and replace with a chromosome 'identifier' that is the chunk nunber
# exclude header 
# make sure to be tab separated! otherwise python script won't know what to do
grep -v "#" $repdir/slim.output.${state}Contraction.${j}.vcf | awk -v chr=$i '{OFS="\t";$1=""; print chr,$0}' >>  $outVCF

done
done
done
done
done