########### THIS SCRIPT (SHELL) WILL MAKE 1 makefile for every sample to be used in paleomix ########
# note you have to run this with source not sh because it uses cd (or submit as a job)
# the Nextseq produces 8 fastq files because there are four 'lanes' (fluidically connected)
# so this accounts for that (Baja samples + LANGEDOG)
############  directories ########
gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
# scripts:
scriptDir=$gitDir/scripts/scripts_20180521/data_processing


########### file locations: ########
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures
fastqs=$wd/fastqs
makefileDir=$scriptDir/paleomixPipeline/makefiles
modernTemplate=$makefileDir/makefile_template.modernDNA.Nextseq.multipleLanes.yaml
modHeaders=$wd/samples/bajaCaptures.plusDogs.txt


mkdir -p $makefileDir/modernMakefiles


########### modern makefiles ########
# modern dna
cat $modHeaders | while read header
do
echo $header
# note the need for double quotation marks for sed
# make a new version of the makefile
/bin/cp $modernTemplate $makefileDir/modernMakefiles/${header}.paleomix.makefile.yaml 
newMake=$makefileDir/modernMakefiles/${header}.paleomix.makefile.yaml
# for now NAME OF TARGET and SAMPLE are going to be the same
sed -i'' "s/NAME_OF_TARGET:/$header:/g" $newMake
sed -i'' "s/NAME_OF_SAMPLE:/$header:/g" $newMake
sed -i'' "s/NAME_OF_LIBRARY:/${header}_1a:/g" $newMake
# for nextseq there are four 'lanes' that were fluidically connected and each generate a fastq
sed -i'' "s/NAME_OF_LANE1:/Lane_1:/g" $newMake
# use different delims (|) to avoid filepath slash confusion:
sed -i'' 's|: PATH_WITH_WILDCARDS1|: '${fastqs}\/${header}_S*_L001_R{Pair}_*fastq.gz'|g' $newMake
sed -i'' "s/NAME_OF_LANE1:/Lane_2:/g" $newMake
# use different delims (|) to avoid filepath slash confusion:
sed -i'' 's|: PATH_WITH_WILDCARDS2|: '${fastqs}\/${header}_S*_L002_R{Pair}_*fastq.gz'|g' $newMake
sed -i'' "s/NAME_OF_LANE1:/Lane_3:/g" $newMake
# use different delims (|) to avoid filepath slash confusion:
sed -i'' 's|: PATH_WITH_WILDCARDS3|: '${fastqs}\/${header}_S*_L003_R{Pair}_*fastq.gz'|g' $newMake
sed -i'' "s/NAME_OF_LANE1:/Lane_4:/g" $newMake
# use different delims (|) to avoid filepath slash confusion:
sed -i'' 's|: PATH_WITH_WILDCARDS4|: '${fastqs}\/${header}_S*_L004_R{Pair}_*fastq.gz'|g' $newMake

# need to say "{Pair}" to have it find both copies. See if this works.
# clear variables
# note that if sed is run on an empty file, it creates randomly named weird empty file. Not a big deal. e.g. sedVeLdud
newMake=''
done
