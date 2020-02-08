# wrapper
# runs fast, can do on home computer
gitdir=/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/
scriptdir=$gitdir/scripts/scripts_20180521/analyses/dadi_inference/grid.search
model=1D.Trim22MSMCModel.PlusExtraEpoch # CA trim22 MSMC model 
script=grid.Search.${model}.dadi.dadiUnits.py # new! 20190911 -- have input be in terms of dadi units and don't rescale
# input directory
genotypeDate=20181119
indir=/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/$genotypeDate/easySFS_projection/neutral/projection-20181221-hetFilter-0.75/dadi-plusMonomorphic

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
# input these from dadi MLE results: CALIFORNIA SPECIFIC PARAMS:
#Nanc=3585
# pop specific out-dir
outdir=/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/$genotypeDate/MSMCModels/$pop/grid.search
mkdir -p $outdir
sfs=$pop-[0-9]*.plusMonomorphic.sfs # the pop specific SFS, skipping [0-9] sample size so you don't have to specify for each pop

#  note that grid is log10 scaled, so will have more data points at the lower values, fewer as you go higher
python $scriptdir/$script --sfs $indir/$sfs --pop $pop  --numGridPoints 100 --nu_Low ${nu_Low} --nu_High ${nu_High} --T_Low ${T_Low} --T_High ${T_High} --outdir $outdir

# gzip output
gzip -f $outdir/dadi.grid.search.$pop.$model.LL.output.txt

# empty variables just in case:
sfs=""
#Nanc=""

################################# then do again with the simplified model #############


# wrapper
# runs fast, can do on home computer
gitdir=/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/
scriptdir=$gitdir/scripts/scripts_20180521/analyses/dadi_inference/grid.search
model=1D.SIMPLIFIED.Trim22MSMCModel.PlusExtraEpoch
 # CA trim22 MSMC model 
script=grid.Search.${model}.dadi.dadiUnits.py # new! 20190911 -- have input be in terms of dadi units and don't rescale
# input directory
genotypeDate=20181119
indir=/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/$genotypeDate/easySFS_projection/neutral/projection-20181221-hetFilter-0.75/dadi-plusMonomorphic

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
# input these from dadi MLE results: CALIFORNIA SPECIFIC PARAMS:
#Nanc=3585
# pop specific out-dir
outdir=/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/$genotypeDate/MSMCModels/$pop/grid.search
mkdir -p $outdir
sfs=$pop-[0-9]*.plusMonomorphic.sfs # the pop specific SFS, skipping [0-9] sample size so you don't have to specify for each pop

#  note that grid is log10 scaled, so will have more data points at the lower values, fewer as you go higher
python $scriptdir/$script --sfs $indir/$sfs --pop $pop  --numGridPoints 100 --nu_Low ${nu_Low} --nu_High ${nu_High} --T_Low ${T_Low} --T_High ${T_High} --outdir $outdir

# gzip output
gzip -f $outdir/dadi.grid.search.$pop.$model.LL.output.txt

# empty variables just in case:
sfs=""
#Nanc=""


