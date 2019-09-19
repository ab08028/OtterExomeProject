############## try threepop and fourpop from treemix ###########

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
wd=/u/flashscratch/a/ab08028/captures/analyses/TREEMIX/$genotypeDate/snp7/fStats # update if using snp7 8 etc 
mkdir -p $wd

######### has ferret #############
k=500
marker="sepCA-BAJ.exclRelatives-testADDINGMFUR-refA2"
infile=${header}.${marker}.frq.strat.treemixFormat.InclFERRET.gz # this has CA-BAJ combined 
#root='CA,BAJ'# root is both
#outdir="root.${root}.mig.${m}.k.${k}.global.${marker}.treemix"
#mkdir -p $wd/$outdir
#cd $outdir
# f3
outdir=/u/flashscratch/a/ab08028/captures/analyses/TREEMIX/20181119/snp7/fStats

threepop -i $treeFileDir/$infile -k ${k} > ${header}.${marker}.f3.out 
# just pull out results with a header
echo -e "pops f3 f3_se zScore" > ${header}.${marker}.f3.out.forR.txt
grep ";" ${header}.${marker}.f3.out  >> ${header}.${marker}.f3.out.forR.txt
# f4
fourpop -i $treeFileDir/$infile -k ${k} > ${header}.${marker}.f4.out 


######## try without ferret: ########
k=500
marker="sepCA-BAJ.exclRelatives"
infile=${header}.${marker}.frq.strat.treemixFormat.gz 
outdir=/u/flashscratch/a/ab08028/captures/analyses/TREEMIX/20181119/snp7/fStats
threepop -i $treeFileDir/$infile -k ${k} > ${header}.${marker}.f3.out 
# just pull out results with a header
echo -e "pops f3 f3_se f3_zScore" > ${header}.${marker}.f3.out.forR.txt
grep ";" ${header}.${marker}.f3.out  >> ${header}.${marker}.f3.out.forR.txt
# f4
fourpop -i $treeFileDir/$infile -k ${k} > ${header}.${marker}.f4.out 
echo -e "pops f4 f4_se f4_zScore" > ${header}.${marker}.f4.out.forR.txt
grep ";" ${header}.${marker}.f4.out  >> ${header}.${marker}.f4.out.forR.txt
