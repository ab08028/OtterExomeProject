
######## concatenate VCF files per replicate, and add a 'chromosome' identifier (that's really a portion identifier)
# loop over all replicates
for i in {1..11}
do
echo "starting replicate $i"
gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptdir=$gitdir/scripts/scripts_20180521/analyses/

model=1D.2Epoch
rundate=20190125 # date of running slim
outdir=/u/flashscratch/a/ab08028/captures/analyses/slim/neutralSimulations/$model/$rundate/replicate_${i} # eventually loop over all replicates 
mkdir -p $outdir/SFS
mkdir /u/flashscratch/a/ab08028/captures/analyses/slim/neutralSimulations/$model/$rundate/allSFSes/ # where SFSes will end up for download to computer
outfile=${model}.rep.${i}.concatted.slim.output.ALL.vcf

# get header from first vcf

echo "concatenating chunks"

grep "#" $outdir/slim.output.1.vcf > $outdir/$outfile
# loop over all chunks
for j in {1..60}
do
# want to select everything but first column (setting it to "") and replace with a chromosome 'identifier' that is the chunk nunber
# exclude header 
# make sure to be tab separated! otherwise python script won't know what to do
grep -v "#" $outdir/slim.output.${j}.vcf | awk -v chr=$i '{OFS="\t";$1=""; print chr,$0}' >> $outdir/$outfile
done


echo "gzipping output"
# gzip the result
gzip -f $outdir/$outfile
# this will give you once vcf for the whole run

# then want to get the SFS from it

# use my SFS script
# this will generate an UNFOLDED SFS in dadi and R format
echo "generating SFS"

python $scriptdir/generate_sfs/make_sfs_without_easySFS/generate1DSFS.py \
--vcf $outdir/${outfile}.gz \
--pop generic \
--outdir $outdir/SFS \
--outPREFIX rep.$j.$model.slim.output


# can then fold it when you do dadi inference.

# and want to gather up all the SFSes
echo "gathering sfses together for download"
/bin/cp -r $outdir/SFS/* /u/flashscratch/a/ab08028/captures/analyses/slim/neutralSimulations/$model/$rundate/allSFSes/
# use bin/cp to overrride cp -i behavior


done


# can then fold it when you do dadi inference.

