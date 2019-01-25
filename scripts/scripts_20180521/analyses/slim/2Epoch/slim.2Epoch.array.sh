#! /bin/bash
#$ -cwd
#$ -l h_rt=02:00:00,h_data=2G
#$ -o /u/flashscratch/a/ab08028/captures/reports/slim
#$ -e /u/flashscratch/a/ab08028/captures/reports/slim
#$ -m abe
#$ -M ab08028
#$ -t 1-60


######### 2 Epoch script generates 100 x 1kb independent blocks ###########
# want a total of 6000 blocks. so 60 instances of this script for one replicate.

rep=$1 # doing one replicate, then will set from command line from submission script
rundate=$2 # date arrays are submitted; set in submitter so as not to have jobs on different days
outdir=$SCRATCH/captures/analyses/slim/$model/$rundate/replicate_${rep} # set this in submission script 
mkdir -p $outdir
######### programs #########
# load the proper gcc 
source /u/local/Modules/default/init/modules.sh
module load gcc/6.3.0
slim=/u/home/a/ab08028/klohmueldata/annabel_data/bin/SLiM/slim_build/slim # location of slim

############## files and dirs ############
gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/ # project github location
scriptdir=$gitdir/scripts/scripts_20180521/analyses/slim/ # location of slim scripts
slimscript=generic.2Epoch.100kb.10genContraction.20180125.slim # specific slim script
model=2Epoch

######## parameters #############
seed=1 # figure out seed ; dont use SGE TASK ID
mu=8.64e-9
r=1e-8
ss=7 # sample size in individuals
nanc=4000 # ancestral size from dadi 
nu=30 # bottleneck size from dadi

#t=10 # time before present that contraction occured
#burnin=$((nanc*10)) # burn in time (nanc *10)
#toutput=$((burnin+contractdur)) # burnin + generations bp that contraction occurs; simulation ends at toutput

# vars: v_NANC  = 4000; v_NU = 30; v_ 10 ; v_MU = 8.64e-9 ; v_SS (sample size) ; v_OUTFILE (path to outfile destination)
# note for sample size :" "A sample of individuals (not genomes, note Ð unlike the outputSample() and outputMSSample() methods) of size sampleSize from the subpopulation will be output.  The sample may be done either with or without replacement, as specified by replace;"
# math calculation. burn in is 10* Nanc, then add contraction time bp duration 
# for now can't get times to work 

$slim \
-long \
-seed $seed \
-d v_MU=$mu \
-d v_R=$r \
-d v_SS=$ss \
-d v_NANC=$nanc \
-d v_NU=$nu \
-d v_OUTFILE="'$outdir'" \
-d block=${SGE_TASK_ID} \
$scriptdir/$slimscript
# script comes at end
# long prints out what vars are set as 
# this will output vcf and full population state.


sleep 10m


####### for testing on home computer ########
# outdir="'/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/analyses/slim/sandbox/'"
# slim \
# -long \
# -seed $seed \
# -d v_MU=$mu \
# -d v_R=$r \
# -d v_SS=$ss \
# -d v_NANC=$nanc \
# -d v_NU=$nu \
# -d v_OUTFILE=$outdir \
# -d v_TCONTRACT=$t \
# -d block=1 \
# generic.2Epoch.100kb.slim
