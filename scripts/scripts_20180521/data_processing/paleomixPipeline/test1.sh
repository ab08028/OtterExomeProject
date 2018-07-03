#TEST:::
# location of github:  which may be on remote server or local laptop
gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
# scripts:
scriptDir=$gitDir/scripts/scripts_20180521/data_processing/paleomix
# script to run: 
scriptname=plmx_Step_1_runPlmx.sh

# file locations:
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures
fastqs=$wd/fastqs
makefileDir=$gitDir/paleomixPipeline/makefiles
# makefile:

# outdirectory:
outdir=$wd/paleomix
mkdir -p $outdir
# job info: 
errorLocation=/u/flashscratch/a/ab08028/captures/reports/paleomix # report location
user=ab08028 # where emails are sent
header=46_Elut_BER_29
qsub -e $errorLocation -o $errorLocation -M $user -N plmx${c} -pe shared 8 \
$scriptDir/$scriptname $makefileDir/${header}.paleomix.makefile.yam $outdir