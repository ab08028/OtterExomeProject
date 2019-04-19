########### THIS SCRIPT (SHELL) WILL MAKE 1 makefile for every sample to be used in paleomix ########
# note you have to run this with source not sh because it uses cd (or submit as a job)
############  directories ########
gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
# scripts:
scriptDir=$gitDir/scripts/scripts_20180521/data_processing


########### file locations: ########
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures
fastqs=$wd/fastqs
makefileDir=$scriptDir/paleomixPipeline/makefiles
#modernTemplate=$makefileDir/makefile_template.modernDNA.yaml
ancientTemplate=$makefileDir/makefile_template.aDNA.fullProcessing.IndelR.ElutMfur.yaml
#ancHeaders=$wd/samples/bestAncientHeaders.txt # 3 best aDNAs 
#modHeaders=$wd/samples/modernSamples.txt # choose 3 CA and 3 AK to process; want to downsample at some point along the way (maybe at BAM stage)
# A13 has _L_ , A29-30 don't


mkdir -p $makefileDir/ancientMakefiles
#mkdir -p $makefileDir/modernMakefiles
########### ancient makefiles ########
# instead of START/END, using while read list of samples
#cat $ancHeaders | while read header
#do
#echo $header
# note the need for double quotation marks for sed
# make a new version of the makefile
######## *** ALWAYS MAKE SURE SETTINGS ARE APPROPRIATE *** ###########
# calling /bin/cp because my cp is aliased to be interactive

for header in A30_Elut_CA_SM_35_SN1_CAP A29_Elut_CA_SM_30_SN2_CAP
do
/bin/cp $ancientTemplate $makefileDir/ancientMakefiles-fullProcessing-usethis/${header}.paleomix.makefile.yaml
newMake=$makefileDir/ancientMakefiles-fullProcessing-usethis/${header}.paleomix.makefile.yaml
# for now NAME OF TARGET and SAMPLE are going to be the same
sed -i'' "s/NAME_OF_TARGET:/$header:/g" $newMake
sed -i'' "s/NAME_OF_SAMPLE:/$header:/g" $newMake
sed -i'' "s/NAME_OF_LIBRARY:/${header}_1a:/g" $newMake
# for now just naming Lane: Lane 1 because just one lane of novaseq
sed -i'' "s/NAME_OF_LANE:/Lane_1:/g" $newMake
# use different delims (|) to avoid filepath slash confusion:
sed -i'' 's|: PATH_WITH_WILDCARDS|: '${fastqs}\/${header}_S*_R{Pair}_*fastq.gz'|g' $newMake

# clear variables
newMake=''
done

for header in A13_Elut_CA_AN_388_SN1_2CAP_screen
do
/bin/cp $ancientTemplate $makefileDir/ancientMakefiles-fullProcessing-usethis/${header}.paleomix.makefile.yaml
newMake=$makefileDir/ancientMakefiles-fullProcessing-usethis/${header}.paleomix.makefile.yaml
# for now NAME OF TARGET and SAMPLE are going to be the same
sed -i'' "s/NAME_OF_TARGET:/$header:/g" $newMake
sed -i'' "s/NAME_OF_SAMPLE:/$header:/g" $newMake
sed -i'' "s/NAME_OF_LIBRARY:/${header}_1a:/g" $newMake
# for now just naming Lane: Lane 1 because just one lane of novaseq
sed -i'' "s/NAME_OF_LANE:/Lane_1:/g" $newMake
# use different delims (|) to avoid filepath slash confusion:
sed -i'' 's|: PATH_WITH_WILDCARDS|: '${fastqs}\/${header}_S*_L*_R{Pair}_*fastq.gz'|g' $newMake

done
########### modern makefiles ########
# modern dna
#cat $modHeaders | while read header
#do
#echo $header
# note the need for double quotation marks for sed
# make a new version of the makefile
#/bin/cp $modernTemplate $makefileDir/modernMakefiles/${header}.paleomix.makefile.yaml
#newMake=$makefileDir/modernMakefiles/${header}.paleomix.makefile.yaml
# for now NAME OF TARGET and SAMPLE are going to be the same
#sed -i'' "s/NAME_OF_TARGET:/$header:/g" $newMake
#sed -i'' "s/NAME_OF_SAMPLE:/$header:/g" $newMake
#sed -i'' "s/NAME_OF_LIBRARY:/${header}_1a:/g" $newMake
# for now just naming Lane: Lane 1 because just one lane of novaseq
#sed -i'' "s/NAME_OF_LANE:/Lane_1:/g" $newMake
# use different delims (|) to avoid filepath slash confusion:
#sed -i'' 's|: PATH_WITH_WILDCARDS|: '${fastqs}\/${header}_S*_R{Pair}_*fastq.gz'|g' $newMake
# need to say "{Pair}" to have it find both copies. See if this works.
# clear variables
# note that if sed is run on an empty file, it creates randomly named weird empty file. Not a big deal. e.g. sedVeLdud
#newMake=''
#done
