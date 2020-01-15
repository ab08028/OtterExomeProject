#! /bin/bash
#$ -cwd
#$ -l h_rt=3:00:00,h_data=2G
#$ -N msmc
#$ -pe shared 16
#$ -o /u/flashscratch/a/ab08028/captures/reports/msmc
#$ -e /u/flashscratch/a/ab08028/captures/reports/msmc
#$ -m abe
#$ -M ab08028


# run MSMC with default settings (for now)

source /u/local/Modules/default/init/modules.sh
module load anaconda
source activate MaCsSimulations # a conda env with python 3.6 in it

msmc_tools=/u/home/a/ab08028/klohmueldata/annabel_data/bin/msmc-tools
msmc=/u/home/a/ab08028/klohmueldata/annabel_data/bin/msmc
rundate=`date +%Y%m%d` # msmc rundate 
OUTDIR=/u/flashscratch/a/ab08028/captures/runMSMCOnSimulations/${rundate}
mkdir -p $OUTDIR

models="dadiModel1.OldShallowContraction  dadiModel2.RecentExtremeContraction"

for model in $models
do
for rep in `seq 1 10`
do
## only once: gather up MSMC files into one dir
INPUTDIR=/u/flashscratch/a/ab08028/captures/analyses/simulateForMSMC/$model/rep_$rep/allMSMCInputFiles
mkdir -p $INPUTDIR
OUTDIR=/u/flashscratch/a/ab08028/captures/runMSMCOnSimulations/$model/rep_$rep/
mkdir -p $OUTDIR
# cp /u/flashscratch/a/ab08028/captures/analyses/simulateForMSMC/$model/rep_$rep/*/*MSMCFormat* $INPUTDIR ## ONLY COPY THINGS OVER ONCE 
# unless you redo the simulation
$msmc -t 16 -o $OUTDIR/msmc.RunOn.${model}.sims.out $INPUTDIR/group_*_block_*.${model}.MSMCFormat.OutputFile.20200114.txt

########## run msmc on the simulated chunks for each replicate #######
done
done
