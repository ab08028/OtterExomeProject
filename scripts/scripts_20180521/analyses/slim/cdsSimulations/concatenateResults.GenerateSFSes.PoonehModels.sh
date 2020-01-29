########### this script will gather up simulation replicates within a single rundate (run again for different rundates)
# for the populations and models you specify
# you can also modify for differnent run dates by supplying the popsModelsRundates list  manually 
# it will gather up summaries of the simulations (one line per mutation) and will concanenate VCF files and make SFSes from them
##### SOMETHING TO LOOK INTO: Jacqueline found that the VCF files from SLIM were kind of weird, which is why she does summaries
# Make sure they are concordant!

# gather up replicates:
numReps=25 # total number of reps you ran 
DesiredReps=20 # how many you'll take out of the 25 (some randomly fail, so picking 20 / 25 for all)
numChunks=20 #
numStates=4 # there are more states now 
checkNumber=$((numChunks*numStates)) # should be 80 as the check number now
# still just going to concat over two states 

# process output of slim 
gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
# choose the specific combination of populations/models/rundates you want? this is awkward... what is best way to do it?
# com is 3epoch (differnt model) 

#rundate=20190424 # set of simulations you're interested in (if is across different rundates you can list popsModelsRundates explicitly)
hs="0 0.5" # set of hs you're interested in
#popMods="AL/1D.2Epoch.1.5Mb.cds AK/1D.2Epoch.1.5Mb.cds" # population and corresponding models you're interested in

# you can have multiple models per population just list as: AK/1D.2Epoch.1.5Mb.cds AK/OtherModel in the popMods variable
scriptdir=$gitdir/scripts/scripts_20180521/analyses/slim/cdsSimulations/
wd=$SCRATCH/captures/analyses/slim/poonehSimulations/ # location of Pooneh's simulations
####### First just do alaska: 
popsModelsRundates="AK/1D.3Epoch.LongerRecovery/20191202/h_0/ AK/1D.3Epoch.LongerRecovery/20191202/h_0.5/"


# for i in $popMods
# do
# for h in $hs
# do
# pm=$i/$rundate/h_$h
# popsModelsRundates=`echo $popsModelsRundates $pm`
# done
# done # this sets up your list of pops, models, rundates and Hs in a list that looks like:
# AK/1D.2Epoch.1.5Mb.cds/20190423/h_0.5/ AL/1D.2Epoch.1.5Mb.cds/20190423/h_0.5/ CA/1D.2Epoch.1.5Mb.cds/20190423/h_0.5/ KUR/1D.2Epoch.1.5Mb.cds/20190423/h_0.5/

# or you can set it up manually if you are using an odd mixture of models/dates --> 
# popsModelsRundates='AK/1D.2Epoch.1.5Mb.cds/20190423/h_0.5/ AL/1D.2Epoch.1.5Mb.cds/20190423/h_0.5/ CA/1D.2Epoch.1.5Mb.cds/20190423/h_0.5/ KUR/1D.2Epoch.1.5Mb.cds/20190423/h_0.5/' # maybe? -- this is kind of awkward, maybe have to deal with diff populations differently?
# not ready yet: COM/1D.3Epoch.1.5Mb.cds/20190423/h_0.5/
# loop through models, populations and 25 replicates


################### first get lists of which runs worked #################

for popsModelsRundate in $popsModelsRundates
do
echo $popsModelsRundate
#pull out the population name for when you make sfses (note that if you change how you label popmodelsrundates you'll have to change how to do this)
pop=${popsModelsRundate%/*/*/*/*} #
i=${popsModelsRundate%/*/*/*} # just the pop/model
############### get top passing runs ############################
# use a little script I made to detect the ones that passed and make lists
# usage script.sh [popModel] [total reps run] [Desired replicates to pull out] [numChunks per replicate to check for completeness] [number of states per replicate (e.g. pre and post contraction)]
chmod +x $scriptdir/FigureOutWhichRunsPassed.sh # make sure it's executable
$scriptdir/FigureOutWhichRunsPassed.poonehSimulations.sh $popsModelsRundate $numReps $DesiredReps $numChunks $numStates # this will make files of the ones that passed 
passingFile=passingReps.FIRST.${DesiredReps}.usethis.txt 
# make sure the file exists and is $DesiredReps long 
if [ ! -f $wd/$popsModelsRundate/$passingFile ]
then
break
fi

##################################################################
# make the slim script from the maker script:
vcfOutDir=$wd/concattedVCFs/$popsModelsRundate/
#summaryOutDir=$wd/concattedSummaries/$popsModelsRundate/
sfsDir=$wd/SFSes/$popsModelsRundate/
sfsScriptDir=$gitdir/scripts/scripts_20180521/analyses/generate_sfs/make_sfs_without_easySFS

mkdir -p $vcfOutDir
#mkdir -p $summaryOutDir
mkdir -p $sfsDir

for state in Pre Post # pre and post contraction
do
cat $wd/$popsModelsRundate/$passingFile | while read i # instead of a for loop, going to read through my passing file. i is replicate_x not just a number 
do 
echo "$i ${state}Contraction"
repdir=$wd/$popsModelsRundate/${i}

#outSummary=$summaryOutDir/${i}.slim.output.${state}Contraction.allConcatted.summary.txt # skip summaries (Pooneh already did)

# vcf with all mutations together:
outVCF=$vcfOutDir/${i}.slim.output.${state}Contraction.allConcatted.vcf

# mutation subtypes vcfs:
mut1VCF=$vcfOutDir/${i}.slim.output.${state}Contraction.allConcatted.mutType1.vcf
# mutation subtypes vcfs:
mut2VCF=$vcfOutDir/${i}.slim.output.${state}Contraction.allConcatted.mutType2.vcf

#vcf headers:
grep "#" $repdir/slim.output.${state}Contraction.1.vcf > $outVCF
grep "#" $repdir/slim.output.${state}Contraction.1.vcf > $mut1VCF
grep "#" $repdir/slim.output.${state}Contraction.1.vcf > $mut2VCF
# summary header:
#grep "replicate" $repdir/slim.output.${state}Contraction.1.summary.txt > $outSummary

# then loop over all chunks
for j in $(seq 1 $numChunks)
do
# concat summaries and gzip
#echo "concatenating summaries"
#  grep -v "^$" removes extra blank lines
#grep -v "replicate" $repdir/slim.output.${state}Contraction.${j}.summary.txt | grep -v "^$" >> $outSummary
# concat vcfs and gzip
#echo "concatenating vcf chunks"
# # want to select everything but first column (setting it to "") and replace with a chromosome 'identifier' that is the chunk nunber
# exclude header 
# make sure to be tab separated! otherwise python script won't know what to do
grep -v "#" $repdir/slim.output.${state}Contraction.${j}.vcf | awk -v chr=$j '{OFS="\t";$1=""; print chr,$0}' >>  $outVCF
# # separate mutation types
# mutation type 1: 
grep -v "#" $repdir/slim.output.${state}Contraction.${j}.vcf | grep "MT=1"  | awk -v chr=$j '{OFS="\t";$1=""; print chr,$0}' >> $mut1VCF
# mutation type 2:
grep -v "#" $repdir/slim.output.${state}Contraction.${j}.vcf | grep "MT=2"  | awk -v chr=$j '{OFS="\t";$1=""; print chr,$0}' >> $mut2VCF

done
# gzip the outputs:
echo "Gzipping results"
#gzip -f $outSummary
gzip -f $outVCF
gzip -f $mut1VCF
gzip -f $mut2VCF
# make SFSes: mut type 1
# to just get population name, not whole model date etc do
echo "making SFSes"
python $sfsScriptDir/generate1DSFS.py \
--vcf ${mut1VCF}.gz \
--pop ${pop}_sim \
--outdir $sfsDir \
--outPREFIX $i.${state}Contraction.slim.mutType1

# make SFSes: mut type 2

python $sfsScriptDir/generate1DSFS.py \
--vcf ${mut2VCF}.gz \
--pop ${pop}_sim \
--outdir $sfsDir \
--outPREFIX $i.${state}Contraction.slim.mutType2
done
done
done

########################## California needs a bit of work: VCFs have 7 individuals sampled (should be 6) #########
# gather up replicates:
numReps=25 # total number of reps you ran 
DesiredReps=20 # how many you'll take out of the 25 (some randomly fail, so picking 20 / 25 for all)
numChunks=20 #
numStates=4 # there are more states now 
checkNumber=$((numChunks*numStates)) # should be 80 as the check number now
# still just going to concat over two states 

# process output of slim 
gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
# choose the specific combination of populations/models/rundates you want? this is awkward... what is best way to do it?
# com is 3epoch (differnt model) 

#rundate=20190424 # set of simulations you're interested in (if is across different rundates you can list popsModelsRundates explicitly)
hs="0 0.5" # set of hs you're interested in
#popMods="AL/1D.2Epoch.1.5Mb.cds AK/1D.2Epoch.1.5Mb.cds" # population and corresponding models you're interested in

# you can have multiple models per population just list as: AK/1D.2Epoch.1.5Mb.cds AK/OtherModel in the popMods variable
scriptdir=$gitdir/scripts/scripts_20180521/analyses/slim/cdsSimulations/
wd=$SCRATCH/captures/analyses/slim/poonehSimulations/ # location of Pooneh's simulations

popsModelsRundates="CA/1D.3Epoch.LongerRecovery/20191013/h_0/ CA/1D.3Epoch.LongerRecovery/20191013/h_0.5/"
for popsModelsRundate in $popsModelsRundates
do
echo $popsModelsRundate
#pull out the population name for when you make sfses (note that if you change how you label popmodelsrundates you'll have to change how to do this)
pop=${popsModelsRundate%/*/*/*/*} #
i=${popsModelsRundate%/*/*/*} # just the pop/model
############### get top passing runs ############################
# use a little script I made to detect the ones that passed and make lists
# usage script.sh [popModel] [total reps run] [Desired replicates to pull out] [numChunks per replicate to check for completeness] [number of states per replicate (e.g. pre and post contraction)]
chmod +x $scriptdir/FigureOutWhichRunsPassed.sh # make sure it's executable
$scriptdir/FigureOutWhichRunsPassed.poonehSimulations.sh $popsModelsRundate $numReps $DesiredReps $numChunks $numStates # this will make files of the ones that passed 
passingFile=passingReps.FIRST.${DesiredReps}.usethis.txt 
# make sure the file exists and is $DesiredReps long 
if [ ! -f $wd/$popsModelsRundate/$passingFile ]
then
break
fi

##################################################################
# make the slim script from the maker script:
vcfOutDir=$wd/concattedVCFs/$popsModelsRundate/
#summaryOutDir=$wd/concattedSummaries/$popsModelsRundate/
sfsDir=$wd/SFSes/$popsModelsRundate/
sfsScriptDir=$gitdir/scripts/scripts_20180521/analyses/generate_sfs/make_sfs_without_easySFS

mkdir -p $vcfOutDir
#mkdir -p $summaryOutDir
mkdir -p $sfsDir

for state in Pre Post # pre and post contraction
do
cat $wd/$popsModelsRundate/$passingFile | while read i # instead of a for loop, going to read through my passing file. i is replicate_x not just a number 
do 
echo "$i ${state}Contraction"
repdir=$wd/$popsModelsRundate/${i}

#outSummary=$summaryOutDir/${i}.slim.output.${state}Contraction.allConcatted.summary.txt # skip summaries (Pooneh already did)

# vcf with all mutations together:
outVCF=$vcfOutDir/${i}.slim.output.${state}Contraction.allConcatted.vcf

# mutation subtypes vcfs:
mut1VCF=$vcfOutDir/${i}.slim.output.${state}Contraction.allConcatted.mutType1.vcf
# mutation subtypes vcfs:
mut2VCF=$vcfOutDir/${i}.slim.output.${state}Contraction.allConcatted.mutType2.vcf

#vcf headers:
grep "##" $repdir/slim.output.${state}Contraction.1.vcf > $outVCF
grep "##" $repdir/slim.output.${state}Contraction.1.vcf > $mut1VCF
grep "##" $repdir/slim.output.${state}Contraction.1.vcf > $mut2VCF

### need to cut off the last column of the final line of the ehader that defines the individuals:
grep -v "##" $repdir/slim.output.${state}Contraction.1.vcf | grep "#" | awk '{$(NF--)=""; print}' | tr -s '\n' >> $outVCF
grep -v "##" $repdir/slim.output.${state}Contraction.1.vcf | grep "#" | awk '{$(NF--)=""; print}' | tr -s '\n' >> $mut1VCF
grep -v "##" $repdir/slim.output.${state}Contraction.1.vcf | grep "#" | awk '{$(NF--)=""; print}' | tr -s '\n' >> $mut2VCF

# summary header:
#grep "replicate" $repdir/slim.output.${state}Contraction.1.summary.txt > $outSummary

# then loop over all chunks
for j in $(seq 1 $numChunks)
do
# concat summaries and gzip
#echo "concatenating summaries"
#  grep -v "^$" removes extra blank lines
#grep -v "replicate" $repdir/slim.output.${state}Contraction.${j}.summary.txt | grep -v "^$" >> $outSummary
# concat vcfs and gzip
#echo "concatenating vcf chunks"
# # want to select everything but first column (setting it to "") and replace with a chromosome 'identifier' that is the chunk nunber
# exclude header 
# make sure to be tab separated! otherwise python script won't know what to do

############ Need to remove i6 (7th individual)#############
############# JUST FOR CALIFORNIA!!!!!!!! Have Pooneh adjust script in future:
# awk '{$(NF--)=""; print}' | tr -s '\n' <-- this subtracts 1 from the number of fields so prints all but last column; tr compresses to make sure there's only one \n at end of line

echo "Note: these concatted VCFs have the 7th individual removed using awk to adjust sample size. Don't trust the original AC annotation.\n SFSes made directly from counts so are okay." > $vcfOutDir/NOTE.ABOUT.THESE.VCFs.Had7thIndRemoved.txt
grep -v "#" $repdir/slim.output.${state}Contraction.${j}.vcf | awk -v chr=$j '{OFS="\t";$1=""; print chr,$0}' | awk '{$(NF--)=""; print}' | tr -s '\n' >>  $outVCF
# # separate mutation types
# mutation type 1: 
grep -v "#" $repdir/slim.output.${state}Contraction.${j}.vcf | grep "MT=1"  | awk -v chr=$j '{OFS="\t";$1=""; print chr,$0}' | awk '{$(NF--)=""; print}' | tr -s '\n' >> $mut1VCF
# mutation type 2:
grep -v "#" $repdir/slim.output.${state}Contraction.${j}.vcf | grep "MT=2"  | awk -v chr=$j '{OFS="\t";$1=""; print chr,$0}' | awk '{$(NF--)=""; print}' | tr -s '\n' >> $mut2VCF

done
# gzip the outputs:
echo "Gzipping results"
#gzip -f $outSummary
gzip -f $outVCF
gzip -f $mut1VCF
gzip -f $mut2VCF
# make SFSes: mut type 1
# to just get population name, not whole model date etc do
echo "making SFSes"
python $sfsScriptDir/generate1DSFS.py \
--vcf ${mut1VCF}.gz \
--pop ${pop}_sim \
--outdir $sfsDir \
--outPREFIX $i.${state}Contraction.slim.mutType1

# make SFSes: mut type 2

python $sfsScriptDir/generate1DSFS.py \
--vcf ${mut2VCF}.gz \
--pop ${pop}_sim \
--outdir $sfsDir \
--outPREFIX $i.${state}Contraction.slim.mutType2
done
done
done