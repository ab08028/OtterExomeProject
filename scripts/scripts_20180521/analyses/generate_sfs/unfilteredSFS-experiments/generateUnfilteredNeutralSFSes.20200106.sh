####### Script to pull out neutral regions from unfiltered VCF file
# and project it
# and see if still get a contraction
### other idea: try unprojected SFSes as well (maybe) -- not that necessary I think. skews in data because of missing data. gutenkunst shows projection better

#### HYPOTHESES: want to still reject the 1 epoch model for unfiltered SFS (though won't care about parameters); and want to hopefully reject
# it for unprojected SFS (we'll see)

########### step 0: choose files to use (choosing snp_2 unfiltered variants)###########
module unload java
module load java/1.8.0_111 # must be 1.8 (if you get "unsupported major minor" you've load 1.7 instead )
module load bedtools
vcfdir=/u/flashscratch/a/ab08028/captures/vcf_filtering/20181119_filtered
snpVCF=snp_2_Filter_TrimAlt_raw_variants.vcf.gz # these are raw snps with NO filters applied (will have a lot of bad calls)
neutralBed=${vcfdir}/bedCoords/all_8_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters.min10kb.fromExon.noCpGIsland.noRepeat.noFish.0based.sorted.merged.useThis.bed
outdir=/u/flashscratch/a/ab08028/captures/analyses/SFS/20181119/unfilteredSFS-experiments

GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta

###### need to extract neutral regions #####

##### step 1: extract neutral regions #####
## extract neutral regions using GATK- based on filtering_Step_1_e-ii_pullOutNeutralVariantsFromVCF.forEasySFS.sh (just doing it here) ##

java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/${snpVCF} \
-o $outdir/neutral.${snpVCF%.gz} \
-L $neutralBed

##### step 2: projection ####### 
## do projection ##

source /u/local/Modules/default/init/modules.sh
module unload python
module load python/2.7
module load R

bgzip=/u/home/a/ab08028/klohmueldata/annabel_data/bin/tabix-0.2.6/bgzip
todaysdate=`date +%Y%m%d`
genotypeDate=20181119
noCallFrac=1.0
popFile=/u/flashscratch/a/ab08028/captures/samples/samplesPop.Headers.forEasySFS.3.20181119.txt
perPopHetFilter=1 # for now only filtering sites that are 100% hets (won't be singletons)


gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/
scriptdir=${gitdir}/analyses/generate_sfs/
easySFS=$scriptdir/easySFS.abModified.3.noInteract.Exclude01Sites.HetFiltering.20181121.py  # this is my modification
# this version of script filters het sites that are excessively heterozygous 
populations="CA,AK,AL,COM,KUR"
projections="12,14,20,34,12" # haploids; updated these values on 20181220

### NOTE: projection values must be in same order as populations are in your popFile (this isn't ideal -- at some point I am going to modify the easySFS script)
# note that order is CA,AK,AL,COM,KUR 

easyoutdir=/u/flashscratch/a/ab08028/captures/analyses/SFS/$genotypeDate/unfilteredSFS-experiments/easySFS/neutral/projection-${todaysdate}-hetFilter-${perPopHetFilter}-unfiltered
mkdir -p $easyoutdir

### NOTE: projection values must be in same order as populations are in your popFile (this isn't ideal -- at some point I am going to modify the easySFS script)
# note that order is CA,AK,AL,COM,KUR 
$easySFS -i $outdir/neutral.${snpVCF%.gz} -p $popFile -a -v --proj $projections -f -o $easyoutdir -maxHetFilter $perPopHetFilter


######### also do it with het filter 0.75 ########
perPopHetFilter=0.75 # for now only filtering sites that are 100% hets (won't be singletons)

easyoutdir=/u/flashscratch/a/ab08028/captures/analyses/SFS/$genotypeDate/unfilteredSFS-experiments/easySFS/neutral/projection-${todaysdate}-hetFilter-${perPopHetFilter}-unfiltered
mkdir -p $easyoutdir

### NOTE: projection values must be in same order as populations are in your popFile (this isn't ideal -- at some point I am going to modify the easySFS script)
# note that order is CA,AK,AL,COM,KUR 
$easySFS -i $outdir/neutral.${snpVCF%.gz} -p $popFile -a -v --proj $projections -f -o $easyoutdir -maxHetFilter $perPopHetFilter

###### NOTE: multi-D SFS takes forever (can manually cancel that part)
######## step 3: dadi inference ######
### run dadi with 1 epoch and 2 epoch 

source /u/local/Modules/default/init/modules.sh
module load python/2.7.13_shared
source /u/home/a/ab08028/env_python2.7.13/bin/activate # activate virtual environment 

SCRATCH=/u/flashscratch/a/ab08028/
gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/ # hoffman
scriptdir=$gitdir/scripts/scripts_20180521/analyses/dadi_inference
mu=8.64411385098638e-09
genotypeDate=20181119 # newer gts
sfsDate=20200106 # projection with 0.75 het filter and these projection values ## for COM only
hetFilters="1 0.75"
todaysdate=`date +%Y%m%d`-unfilteredSFS
captures=$SCRATCH/captures/
dadidir=$captures/analyses/dadi_inference/unfilteredSFS-experiments
mkdir -p $dadidir
sfssuffix=sfs # not using monomorphic bin don't care
scripts='1D.1Epoch.dadi.py 1D.2Epoch.dadi.py'

for hetFilter in $hetFilters
do
sfsdir=$captures/analyses/SFS/$genotypeDate/unfilteredSFS-experiments/easySFS/neutral/projection-${sfsDate}-hetFilter-${hetFilter}-unfiltered/dadi/
for pop in CA AK AL COM KUR
do


# get total sites from total sites file that was written out as part of my easySFS scripts

######### THIS IS A STAND-IN L from old projections -- should be ballpark correct but dont' actually care about exact parameters
# just want scripts to run. but this isnt' the perfect unfiltered L!!! This L is too small!!! But just care about model fit and dadi untis

L=`grep $pop /u/flashscratch/a/ab08028/captures/analyses/SFS/20181119/easySFS/neutral/projection-20181221-hetFilter-0.75/dadi-plusMonomorphic/$pop-[0-9]*.totalSiteCount.L.withMonomorphic.txt | awk '{print $2}'` ## not doing this. What L's to use?
## NOTE THESE ARE STAND IN Ls not REAL Ls for unfiltered SFS!!!! ###

for script in $scripts
do
model=${script%.dadi.py}
echo "starting inference for $pop for model $model"
outdir=$dadidir/$genotypeDate/$pop/inference_${todaysdate}/hetFilter-${hetFilter}/$model/
mkdir -p $outdir
# carry out inference with 50 replicates that start with different p0 perturbed params:
for i in {1..50}
do
echo "carrying out inference $i for model $model for pop $pop" 
# [0-9] indicates that it's a number, but not specific about proj value
python $scriptdir/$script --runNum $i --pop $pop --mu $mu --L $L --sfs ${sfsdir}/${pop}-[0-9]*.${sfssuffix} --outdir $outdir
done


echo "concatenating results" # note: some models output date and some don't! 
grep rundate -m1 $outdir/${pop}.dadi.inference.${model}.runNum.1.*output > $outdir/${pop}.dadi.inference.${model}.all.output.concatted.txt
for i in {1..50}
do
grep rundate -A1 $outdir/${pop}.dadi.inference.${model}.runNum.${i}.*output | tail -n1 >> $outdir/${pop}.dadi.inference.${model}.all.output.concatted.txt
done
done 
done
done

##################################### try with unprojected SFS ; need to fold ############
##### dadi can fold inside. Do I have check in my script for that? I might! aha -- my dadi script actually does a check for folding and folds! 
 # so don't need to pre-fold the sfs ; cool
source /u/local/Modules/default/init/modules.sh
module load python/2.7.13_shared
source /u/home/a/ab08028/env_python2.7.13/bin/activate # activate virtual environment 

SCRATCH=/u/flashscratch/a/ab08028/
gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/ # hoffman
scriptdir=$gitdir/scripts/scripts_20180521/analyses/dadi_inference
mu=8.64411385098638e-09
genotypeDate=20180806 # older gts, for now
todaysdate=`date +%Y%m%d`-unprojectedOldSFS
captures=$SCRATCH/captures/
dadidir=$captures/analyses/dadi_inference/unprojectedSFS-experiments
mkdir -p $dadidir
scripts='1D.2Epoch.dadi.py'

sfsdir=$captures/analyses/SFS/$genotypeDate/neutralSFS/
for pop in CA AK AL COM KUR
do


# get total sites from total sites file that was written out as part of my easySFS scripts

######### THIS IS A STAND-IN L from old projections -- should be ballpark correct but dont' actually care about exact parameters
# just want scripts to run. but this isnt' the perfect unfiltered L!!! This L is too small!!! But just care about model fit and dadi untis

L=`grep $pop /u/flashscratch/a/ab08028/captures/analyses/SFS/20181119/easySFS/neutral/projection-20181221-hetFilter-0.75/dadi-plusMonomorphic/$pop-[0-9]*.totalSiteCount.L.withMonomorphic.txt | awk '{print $2}'` ## not doing this. What L's to use?
## NOTE THESE ARE STAND IN Ls not REAL Ls for unfiltered SFS!!!! ###

for script in $scripts
do
model=${script%.dadi.py}
echo "starting inference for $pop for model $model"
outdir=$dadidir/$genotypeDate/$pop/inference_${todaysdate}/$model/
mkdir -p $outdir
# carry out inference with 50 replicates that start with different p0 perturbed params:
for i in {1..50}
do
echo "carrying out inference $i for model $model for pop $pop" 
### NOTE: this SFS isn't folded, but the dadi script has a check for folded in it, and will fold inside the script. so don't need to fold ahead of time (I will double check that this works though)
python $scriptdir/$script --runNum $i --pop $pop --mu $mu --L $L --sfs ${sfsdir}/${pop}.all_9.unfolded.sfs.dadi.format.20181105.txt --outdir $outdir
done


echo "concatenating results" # note: some models output date and some don't! 
grep rundate -m1 $outdir/${pop}.dadi.inference.${model}.runNum.1.*output > $outdir/${pop}.dadi.inference.${model}.all.output.concatted.txt
for i in {1..50}
do
grep rundate -A1 $outdir/${pop}.dadi.inference.${model}.runNum.${i}.*output | tail -n1 >> $outdir/${pop}.dadi.inference.${model}.all.output.concatted.txt
done
done 
done

