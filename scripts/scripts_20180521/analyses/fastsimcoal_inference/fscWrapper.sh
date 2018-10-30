### fsc wrapper

# need to input sample size in haploids into file (Pooneh can do manually) -- is equal to folded sfs dim -1  *2

# dir structure
# fsc_inference
# 		|-CA
#			|-1D.1Bottleneck 
#				|- .est, .tpl files
#		|-AK
# 		|-AL
#		|-COM
#		|-KUR
############# program ############
fsc=/u/home/a/ab08028/bin/fsc251

############# parameters #############
model=1D.1Bottleneck
pop="CA"  # eventually do for all populations
genotypeDate=20180806 # date genotypes were called
sfsDate=20181019 # date sfses were made
header=${model}.fsc
############ file structure ############
wd=/u/flashscratch/a/ab08028/captures/analyses/ #working dir
infDir=$wd/fastsimcoal_inference
outdir=$infDir/$pop/$model
mkdir -p $outdir
sfsDir=$wd/SFS/$genotypeDate/neutralSFS/

########### get sfs into inference directory and rename to match .est and .tpl files #########
cp $sfsDir/${pop}_sfs_${sfsDate}_MAFpop0.obs $infDir/${header}_MAFpop0.obs
$fsc -t ${header}.tpl -n100000 -N100000 -m -e ${header}.est -M 0.001 Ðl 10 -L 40 Ðc6 -q