
######## concatenate VCF files per replicate, and add a 'chromosome' identifier (that's really a portion identifier)
########### set up dir structure ##########
rundate=20190227 # date of running slim
model=1D.2Epoch.30gen
pop=generic
wd=/u/flashscratch/a/ab08028/captures/analyses/slim/neutralSimulations/$model/$rundate/ # eventually loop over all replicates 
sfsdir=$wd/allSFSes
mkdir -p $sfsdir
# loop over all replicates
for i in {1..11}
do

echo "starting replicate $i"

gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptdir=$gitdir/scripts/scripts_20180521/analyses/

repdir=$wd/replicate_${i}
concatFile=${model}.rep.${i}.concatted.slim.output.ALL.vcf

echo "concatenating vcf chunks"
# get header from first vcf
grep "#" $repdir/slim.output.1.vcf > $repdir/$concatFile

# loop over all chunks
for j in {1..60}
do
# want to select everything but first column (setting it to "") and replace with a chromosome 'identifier' that is the chunk nunber
# exclude header 
# make sure to be tab separated! otherwise python script won't know what to do
grep -v "#" $repdir/slim.output.${j}.vcf | awk -v chr=$i '{OFS="\t";$1=""; print chr,$0}' >> $repdir/$concatFile
done

echo "gzipping output"

# gzip the result
gzip -f $repdir/$concatFile
# this will give you once vcf for the whole run

# then want to get the SFS from it

# use my SFS script
# this will generate an UNFOLDED SFS in dadi and R format

echo "generating SFS"

python $scriptdir/generate_sfs/make_sfs_without_easySFS/generate1DSFS.py \
--vcf $repdir/${concatFile}.gz \
--pop $pop \
--outdir $sfsdir \
--outPREFIX rep.$i.$model.slim.output


# can then fold it when you do dadi inference.

done





