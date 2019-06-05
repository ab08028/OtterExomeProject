######## can run in the shell (is very fast) ########
# after a lot of testing of treemix paramters (see sandbox scripts)
# I have found that the paramters that yield the most sensical tree are:
# -k 500 (or any amount of ld pruning); -global (from Gauntiz paper); use sample size correction (don't disable with -noss); 


#### Also want to use files that have excluded 3 close relatives from the data set to not skew allele frequencies 
#### Also want to do with and without baja 

module load treemix
gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/
scriptdir=$gitdir/analyses/TREEMIX/
genotypeDate=20181119
vcfdir=/u/flashscratch/a/ab08028/captures/vcf_filtering/${genotypeDate}_filtered/
treeFileDir=$vcfdir/treemixFormat/
header=snp_7_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants
### using snp7 because it contains admixed individuals
### and I want those migration edges to show up
### it does also contain relatives which isn't amazing but shouldn't make a huge difference 
### can also do with snp9a and see the difference


#header=snp_9a_forPCAetc_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants
wd=/u/flashscratch/a/ab08028/captures/analyses/TREEMIX/$genotypeDate/snp7 # update if using snp7 8 etc 
mkdir -p $wd

######### With BAJA (no relatives) #############
k=500
marker="sepCA-BAJ.exclRelatives"
infile=${header}.${marker}.frq.strat.treemixFormat.gz # this has CA-BAJ combined 
root='CA'
for m in {0..10}
do
outdir="root.${root}.mig.${m}.k.${k}.global.${marker}.treemix"
mkdir -p $wd/$outdir
treemix -i $treeFileDir/$infile -m ${m} -root ${root} -k ${k} -global -o $wd/$outdir/$outdir
done


######### CA Only (no Baja, no relatives) #############
k=500
marker="noBAJA.exclRelatives"
infile=${header}.${marker}.frq.strat.treemixFormat.gz # this has CA-BAJ separated 
root='CA'
for m in {0..10}
do
outdir="root.${root}.mig.${m}.k.${k}.global.${marker}.treemix"
mkdir -p $wd/$outdir
treemix -i $treeFileDir/$infile -m ${m} -root ${root} -k ${k} -global -o $wd/$outdir/$outdir
done