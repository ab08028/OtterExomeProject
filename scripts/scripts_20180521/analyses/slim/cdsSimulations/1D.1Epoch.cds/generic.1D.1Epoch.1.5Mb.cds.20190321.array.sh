#! /bin/bash
#$ -cwd
#$ -l h_rt=10:00:00,h_data=2G
#$ -o /u/flashscratch/a/ab08028/captures/reports/slim
#$ -e /u/flashscratch/a/ab08028/captures/reports/slim
#$ -m abe
#$ -M ab08028
#$ -t 1-14

# going to simulate 1.5MB X 14 (=21MB)
# each 1.5Mb 'chromosome' will have 1000 genes
# each gene is 1500bp long and has an internal 1e-08 recomb rate
# with 1e-03 recomb rate between genes
# genes can have NS or S mutations with ratio 2.31:1
# S mutations are neutral
# NS mutations are drawn from bernard gamma dist DFE, and can have any dom coefficient (set as a variable in this script)
# this script produces a VCF file which you can parse 
# not sure if I want fixed mutations to be included or not. I kind of think I do want them because I count them in my empirical data

######### 2 Epoch script generates 100 x 1kb independent blocks ###########
# want a total of 6000 blocks. so 60 instances of this script for one replicate.
model=1D.1Epoch.cds
rep=$1 # doing one replicate, then will set from command line from submission script
rundate=$2 # date arrays are submitted; set in submitter so as not to have jobs on different days
wd=$SCRATCH/captures/analyses/slim/cdsSimulations/$model/$rundate/
outdir=$wd/replicate_${rep} # set this in submission script 
mkdir -p $outdir
######### programs #########
# load the proper gcc 
source /u/local/Modules/default/init/modules.sh
module load gcc/6.3.0
slim=/u/home/a/ab08028/klohmueldata/annabel_data/bin/SLiM/slim_build/slim # location of slim

############## files and dirs ############
gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/ # project github location
scriptdir=$gitdir/scripts/scripts_20180521/analyses/slim/cdsSimulations/$model # location of slim scripts
slimscript=generic.1D.1Epoch.1.5Mb.cds.20190321.slim # specific slim script


######## parameters #############
todaysdate=`date +%Y%m%d`

# this seed will account for runs run on differet days, and if they start during the same second there will also be adjustments for the run number/task id
# it is not perfectly random, but in general few runs will set the seed at the exact time and the exactness is pretty good (598928/600000 seeds were unique when tested with a for loop)
# in real life, the array jobs will be starting at lots of different times, so mostly it isn't an issue. just a small percentage might start at exact same time
seed=$(($todaysdate+$RANDOM+(($RANDOM*$rep*10))+$SGE_TASK_ID)) # uses date, plus random , plus the replicate and SGE task id. so no task should be the same, and no jobs run on different days should be the same, even if they have same task id
# so if two tasks with the same task id across different reps start in the same second (get same random), then they will still be different because of 10*rep 
mu=8.64e-9
#r=1e-8 -- setting r inside the script; 1e-08 within gene, 1e-3 between genes 
ss=7 # sample size in individuals
nanc=4000 # ancestral size from dadi 
domH=0 # dominance coefficient for NS mutations
echo "mu:$mu ss:$ss nanc:$nanc" > $wd/params.txt

$slim \
-long \
-seed $seed \
-d v_MU=$mu \
-d v_SS=$ss \
-d v_NANC=$nanc \
-d v_H=$domH \
-d v_OUTFILE="'$outdir'" \
-d chunk=${SGE_TASK_ID} \
$scriptdir/$slimscript
# script comes at end
# long prints out what vars are set as 
# this will output vcf and full population state.


sleep 10m



