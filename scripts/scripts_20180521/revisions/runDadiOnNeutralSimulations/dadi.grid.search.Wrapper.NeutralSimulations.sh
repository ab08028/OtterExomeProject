#! /bin/bash
#$ -cwd
#$ -l h_rt=24:00:00,h_data=8G
#$ -N dadi_grid
#$ -o /u/flashscratch/a/ab08028/captures/reports/dadi
#$ -e /u/flashscratch/a/ab08028/captures/reports/dadi
#$ -m abe
#$ -M ab08028
# wrapper
source /u/local/Modules/default/init/modules.sh
module load python/2.7.13_shared
source /u/home/a/ab08028/env_python2.7.13/bin/activate # activate virtual environment 

#SCRATCH=/u/flashscratch/a/ab08028/
#captures=$SCRATCH/captures
gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/ # hoffman
scriptdir=$gitdir/scripts/scripts_20180521/analyses/dadi_inference/grid.search

model=1D.2Epoch
todaysdate=`date +%Y%m%d`

script=grid.Search.${model}.dadi.dadiUnits.py # new! 20190911 -- have input be in terms of dadi units and don't rescale
# input directory
#genotypeDate=20181119
sfsdir=/u/home/a/ab08028/klohmueldata/annabel_data/captures/analyses/slim/neutralSimulations_containsTarballs/neutralSimsForManuscript_20200511/CA.1D.2Epoch.35Gen.200Inds/20200224/allSFSes
dadidir=/u/home/a/ab08028/klohmueldata/annabel_data/captures/analyses/dadi_inference/

simulationModel='CA.1D.2Epoch.35Gen.200Inds.postContraction'
################ search same points for all populations #############
### as of 20190905 I want to switch this to dadi units instead of "real" units --

# old version with "real" units: 
#nu_Low=0.8 # 
#nu_High=6000 # fastsimcoal was 300 or so, but go way higher 
#T_Low=0.1 # 
#T_High=6000 # fastsimcoal was 70 or so, but go higher
# but Nanc differs for each population

############# dadi units : ############ diving by 5000 to get the dadi units
nu_Low=0.00016  # diving 0.8 by 5000 to get the dadi units
nu_High=1.2 # diving 6000 by 5000 to get the dadi units
T_Low=0.00001 # dividing by 2*5000 to get dadi units
T_High=0.6
################# California #########################
pop=CA
for rep in 1 2 3 4 5 6 7 8 9 11
do
outdir=$dadidir/runningDadiOnNeutralSimulations/$simulationModel/grid_search_$todaysdate/$model/replicate_${rep}
mkdir -p $outdir
sfs=$sfsdir/$pop.replicate_${rep}.${simulationModel}.slim.output.unfolded.sfs.dadi.format.txt # the pop specific SFS, skipping [0-9] sample size so you don't have to specify for each pop

#  note that grid is log10 scaled, so will have more data points at the lower values, fewer as you go higher
python $scriptdir/$script --sfs $sfs --pop $pop  --numGridPoints 100 --nu_Low ${nu_Low} --nu_High ${nu_High} --T_Low ${T_Low} --T_High ${T_High} --outdir $outdir

# gzip output
gzip -f $outdir/dadi.grid.search.$pop.$model.LL.output.txt

done

deactivate # deactivate virtualenv
