#### This is a master script generator for the read mapping pipeline. It fills
# in the locations of the programs you use and your output directories. You can run it and it will
# generate the scripts appropriate to your set up.
# do this differently for modern and ancient
######## 
date=`date +%Y%m%d`
wd=/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/
scriptdir=$wd/scripts/generic/ # where your scripts are
project_specific=$wd/scripts/project-specific/scripts_${date}
mkdir -p $scriptdir
mkdir -p ${project_specific} # where you want new version to be
# copy all scripts to that new location:
cp -r $scriptdir/* $project_specific
######## program locations ############

picardlocation=/u/home/a/ab08028/klohmueldata/annabel_data/bin/Picard_2.8.1/picard.jar


###### directory locations #########
scratchLocation=
workingDirectory=

###### submission fields #########
username=ab08028
reportLocation= # this is where -o and -e will go





# update scripts

for file in `find $project_specific -type f -iname "*.sh"`
do
sed -i.'' "s/picardLocation/$picardLocation/g" scriptdir/*/*.sh 



