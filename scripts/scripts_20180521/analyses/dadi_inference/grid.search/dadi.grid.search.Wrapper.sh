# wrapper
# runs fast, can do on home computer
gitdir=/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/
scriptdir=$gitdir/scripts/scripts_20180521/analyses/dadi_inference/grid.search
model=1D.2Epoch
script=grid.Search.${model}.dadi.py
# input directory
genotypeDate=20181119
indir=/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/$genotypeDate/easySFS_projection/neutral/projection-20181221-hetFilter-0.75/dadi-plusMonomorphic

################ search same points for all populations #############
nu_Low=0.8 # 
nu_High=6000 # fastsimcoal was 300 or so, but go way higher 
T_Low=0.1 # 
T_High=6000 # fastsimcoal was 70 or so, but go higher
# but Nanc differs for each population
# so when I plot, should I rescale things? hm.
################# California #########################
pop=CA
# input these from dadi MLE results: CALIFORNIA SPECIFIC PARAMS:
Nanc=3585
# pop specific out-dir
outdir=/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/$genotypeDate/grid.search/$pop
mkdir -p $outdir
sfs=$pop-[0-9]*.plusMonomorphic.sfs # the pop specific SFS, skipping [0-9] sample size so you don't have to specify for each pop

#  note that grid is log10 scaled, so will have more data points at the lower values, fewer as you go higher
python $scriptdir/$script --sfs $indir/$sfs --pop $pop --Nanc $Nanc --numGridPoints 100 --nu_Low ${nu_Low} --nu_High ${nu_High} --T_Low ${T_Low} --T_High ${T_High} --outdir $outdir

# gzip output
gzip -f $outdir/dadi.grid.search.$pop.$model.LL.output.txt

# empty variables just in case:
sfs=""
Nanc=""
################# Alaska #########################
pop=AK
# input these from dadi MLE results: CALIFORNIA SPECIFIC PARAMS:
Nanc=4381
outdir=/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/$genotypeDate/grid.search/$pop
mkdir -p $outdir
sfs=$pop-[0-9]*.plusMonomorphic.sfs # the pop specific SFS, skipping [0-9] sample size so you don't have to specify for each pop

#  note that grid is log10 scaled, so will have more data points at the lower values, fewer as you go higher
python $scriptdir/$script --sfs $indir/$sfs --pop $pop --Nanc $Nanc --numGridPoints 100 --nu_Low ${nu_Low} --nu_High ${nu_High} --T_Low ${T_Low} --T_High ${T_High} --outdir $outdir

# gzip output
gzip -f $outdir/dadi.grid.search.$pop.$model.LL.output.txt

# empty variables just in case:
sfs=""
Nanc=""


################# Kuril #########################
pop=KUR
# input these from dadi MLE results: CALIFORNIA SPECIFIC PARAMS:
Nanc=4349
outdir=/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/$genotypeDate/grid.search/$pop
mkdir -p $outdir
sfs=$pop-[0-9]*.plusMonomorphic.sfs # the pop specific SFS, skipping [0-9] sample size so you don't have to specify for each pop

#  note that grid is log10 scaled, so will have more data points at the lower values, fewer as you go higher
python $scriptdir/$script --sfs $indir/$sfs --pop $pop --Nanc $Nanc --numGridPoints 100 --nu_Low ${nu_Low} --nu_High ${nu_High} --T_Low ${T_Low} --T_High ${T_High} --outdir $outdir

# gzip output
gzip -f $outdir/dadi.grid.search.$pop.$model.LL.output.txt

# empty variables just in case:
sfs=""
Nanc=""


################# AL #########################
pop=AL
# input these from dadi MLE results: CALIFORNIA SPECIFIC PARAMS:
Nanc=4834
outdir=/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/$genotypeDate/grid.search/$pop
mkdir -p $outdir
sfs=$pop-[0-9]*.plusMonomorphic.sfs # the pop specific SFS, skipping [0-9] sample size so you don't have to specify for each pop

#  note that grid is log10 scaled, so will have more data points at the lower values, fewer as you go higher
python $scriptdir/$script --sfs $indir/$sfs --pop $pop --Nanc $Nanc --numGridPoints 100 --nu_Low ${nu_Low} --nu_High ${nu_High} --T_Low ${T_Low} --T_High ${T_High} --outdir $outdir

# gzip output
gzip -f $outdir/dadi.grid.search.$pop.$model.LL.output.txt

# empty variables just in case:
sfs=""
Nanc=""

############## skipping commanders for now, because it's weird ################
