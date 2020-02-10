############## Want to try ZooROH #############
####
# need to convert to "Oxford Gen" format (pg 10 of manual)
module load plink
### can be run in the shell ()

##### want to use existing ped files that i generated in /Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/analyses/ROH_plink/runROH.plink.sh

pedDir=/u/flashscratch/a/ab08028/captures/analyses/ROH_plink
zooDir=/u/flashscratch/a/ab08028/captures/analyses/zooROH
### to make files easier ot work with, insist on maf = 0.001 so that you aren't stuck with monomorphic; this shouldn't remove any snps, just mono 1/1 sites.
# note: zooroh doesn't consider monomorphic sites anyway, but it makes file sizes much smaller. 
for pop in CA AK AL KUR COM
do
popPedDir=$pedDir/$pop
pedFile=plinkFormat.${pop}
outdir=$zooDir/$pop
mkdir -p $outdir

# convert to oxford gen format:
plink --file $popPedDir/$pedFile \
--recode oxford \
--out $outdir/${pop}.Oxford.maf.001 \
--allow-extra-chr \
--const-fid \
--maf 0.001 # this shouldn't get rid of any snps, just get rid of monomorphic sites (--mac not well tested yet says Plink)
# gzip the file afterward
gzip -f $outdir/${pop}.Oxford.maf.001.gen
# need to remove -nan

# get a sample ID file 
awk '{print $2}' $outdir/${pop}.Oxford.maf.001.nosex > $outdir/${pop}.Oxford.maf.001.SampleIDs
# note: don't need to remove -nans - that is only if you are converting from vcf. Plink has already made those missing genotypes 0 0 0 (triple zero is missing)
done
# this corresponds to "GP" format
# but note that phred scores are either 0 or 1, not in between, because it's coming from vcf. Not sure if this is ok or not. Probably since it's high coverage and filtered
