# process output of slim 
gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
models='1D.2Epoch.1.5Mb.cds'
# choose the specific combination of populations/models/rundates you want? this is awkward... what is best way to do it?
# com is 3epoch (differnt model) 
#popsModelsRundates='COM/1D.3Epoch.1.5Mb.cds/20190404/h_0/ AK/1D.2Epoch.1.5Mb.cds/20190404/h_0/ AL/1D.2Epoch.1.5Mb.cds/20190404/h_0 CA/1D.2Epoch.1.5Mb.cds/20190404/h_0 KUR/1D.2Epoch.1.5Mb.cds/20190404/h_0' # maybe? -- this is kind of awkward, maybe have to deal with diff populations differently?
popsModelsRundates='COM/1D.3Epoch.1.5Mb.cds/20190423/h_0.5/ AK/1D.2Epoch.1.5Mb.cds/20190423/h_0.5/ AL/1D.2Epoch.1.5Mb.cds/20190423/h_0.5/ CA/1D.2Epoch.1.5Mb.cds/20190423/h_0.5/ KUR/1D.2Epoch.1.5Mb.cds/20190423/h_0.5/' # maybe? -- this is kind of awkward, maybe have to deal with diff populations differently?
#populations='AK AL CA COM KUR'
#populations="AK AL CA KUR" 
#populations="AK"
# loop through models, populations and 25 replicates
scriptdir=$gitdir/scripts/scripts_20180521/analyses/slim/cdsSimulations/

rundate=`date +%Y%m%d`

#for pop in $populations
#do
#echo $pop
for popsModelsRundate in $popsModelsRundates
do
echo $popsModelsRundate
# make the slim script from the maker script:
wd=$SCRATCH/captures/analyses/slim/cdsSimulations/ # combine model and rundate
vcfOutDir=$wd/concattedVCFs/$popsModelsRundate
summaryOutDir=$wd/concattedSummaries/$popsModelsRundate
sfsDir=$wd/SFSes/$popsModelsRundate
sfsScriptDir=$gitdir/scripts/scripts_20180521/analyses/generate_sfs/make_sfs_without_easySFS

mkdir -p $vcfOutDir
mkdir -p $summaryOutDir
mkdir -p $sfsDir

for state in Pre Post # pre and post contraction
do
for i in {1..1}
do
echo "rep $i ${state}Contraction"
repdir=$SCRATCH/captures/analyses/slim/cdsSimulations/$popsModelsRundate/replicate_${i}

outSummary=$summaryOutDir/rep.${i}.slim.output.${state}Contraction.allConcatted.summary.txt

# vcf with all mutations together:
outVCF=$vcfOutDir/rep.${i}.slim.output.${state}Contraction.allConcatted.vcf

# mutation subtypes vcfs:
mut1VCF=$vcfOutDir/rep.${i}.slim.output.${state}Contraction.allConcatted.mutType1.vcf
# mutation subtypes vcfs:
mut2VCF=$vcfOutDir/rep.${i}.slim.output.${state}Contraction.allConcatted.mutType2.vcf

#vcf headers:
grep "#" $repdir/slim.output.${state}Contraction.1.vcf > $outVCF
grep "#" $repdir/slim.output.${state}Contraction.1.vcf > $mut1VCF
grep "#" $repdir/slim.output.${state}Contraction.1.vcf > $mut2VCF
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
# # separate mutation types
# mutation type 1: 
grep -v "#" $repdir/slim.output.${state}Contraction.${j}.vcf | grep "MT=1"  | awk -v chr=$i '{OFS="\t";$1=""; print chr,$0}' >> $mut1VCF
# mutation type 2:
grep -v "#" $repdir/slim.output.${state}Contraction.${j}.vcf | grep "MT=2"  | awk -v chr=$i '{OFS="\t";$1=""; print chr,$0}' >> $mut2VCF

done
# gzip the outputs:
echo "Gzipping results"
gzip -f $outSummary
gzip -f $outVCF
gzip -f $mut1VCF
gzip -f $mut2VCF
# make SFSes: mut type 1
# to just get population name, not whole model date etc do
pop=${popsModelsRundate%/*/*}
echo "making SFSes"
python $sfsScriptDir/generate1DSFS.py \
--vcf ${mut1VCF}.gz \
--pop ${pop}_sim \
--outdir $sfsDir \
--outPREFIX rep.$i.${state}Contraction.slim.mutType1

# make SFSes: mut type 2

python $sfsScriptDir/generate1DSFS.py \
--vcf ${mut2VCF}.gz \
--pop ${pop}_sim \
--outdir $sfsDir \
--outPREFIX rep.$i.${state}Contraction.slim.mutType2
done
done
done
