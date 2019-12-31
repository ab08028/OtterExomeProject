# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 16:46:27 2018

@author: annabelbeichman
"""
import matplotlib
matplotlib.use('Agg') # so graphics show up on hoffman
import sys
import argparse
import dadi
from dadi import Numerics, PhiManip, Integration, Spectrum
from numpy import array # don't comment this out
import datetime
todaysdate=datetime.datetime.today().strftime('%Y%m%d')


modelName="1D.1Bottleneck.FixTRec"

############### Parse input arguments ########################
parser = argparse.ArgumentParser(description='Infer a '+ modelName +' model from a 1D folded SFS in dadi')
parser.add_argument("--runNum",required=True,help="iteration number (e.g. 1-50)")
parser.add_argument("--pop",required=True,help="population identifier, e.g. 'CA'")
parser.add_argument("--mu",required=True,help="supply mutation rate in mutation/bp/gen")
parser.add_argument("--L",required=True,help="number of called neutral sites that went into making SFS (monomorphic+polymorphic)")
parser.add_argument("--sfs",required=True,help="path to folded SFS in dadi format that is output from easySFS pipeline (mask optional)")
parser.add_argument("--outdir",required=True,help="path to output directory")
# usage:
# python 1D.Bottleneck.dadi.py --runNum $i --pop CA --mu 8.64411385098638e-09 --L 4193488 --sfs [path to sfs] --outdir [path to outdir]
args = parser.parse_args()
runNum=str(args.runNum)
pop=str(args.pop)
mu=float(args.mu)
L=float(args.L)
outdir=str(args.outdir)
sfs=str(args.sfs)
maxiter=100
############### Input data ####################################
# for testing purposes:
#sfs="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/20180806/neutralSFS/CA.all_9.unfolded.sfs.dadi.format.20181105.txt"
fs=dadi.Spectrum.from_file(sfs) # this is folded from easySFS
# fold the fs:
# fs=fs.fold() # folded
# check if it's folded, if not folded, fold it
if fs.folded==False:
    fs=fs.fold()
else:
    fs=fs
############### Set up General Dadi Parameters ########################
ns = fs.sample_sizes # get sample size from SFS (in haploids)
pts_l = [ns[0]+5,ns[0]+15,ns[0]+25] # this should be slightly larger (+5) than sample size and increase by 10
############### Set up Specific Model -- this will change from script to script ########################
# Bernard says to fix bottleneck duration
# 20 generations ~ 20/(3000*2)
def bottleneck(params, ns, pts):
    nuB,nuF,TB = params
    xx = Numerics.default_grid(pts) # sets up grid
    phi = PhiManip.phi_1D(xx) # sets up initial phi for population
    phi = Integration.one_pop(phi, xx, TB, nuB)  # bottleneck
    phi = Integration.one_pop(phi, xx, 0.001, nuF) # recovery
    fs = Spectrum.from_phi(phi, ns, (xx,))
    return fs
param_names=("nuB","nuF","TB")

# 20181024 changing upper bound on pop sizes to 10 because know it doesn't grow up to 100* Nanc; and lowering lower bounds to 1e-4 ; see what happens
# changing starting position to .1, tb to 0.005; TF to 0.001
# 20181031: note from Bernard:
#1) for upper_bound, i'd shrink the upper bound of timesby a lot
#not sure about otters, but for humans i think i used like 0.1
#in any case, even T=1 seems like a really long time for these upper bounds assuming the otter #popn size changes you're modeling are on anthropogenic timescales. i don't know how this works #out given otter generation times and the time scale of their decline though (edited)
#another reason they tend to be smaller is that these are times for that period rather than #times from the present
#2) the times for the lower bounds are supposed to be non-zero
#IIRC if you give a min of 0 things get weird because of collapsing epochs
#i'd throw something like 1e-5
upper_bound = [10, 2,2]
lower_bound = [1e-4, 1e-4,1e-4]
p0 = [0.01,0.1,0.005] # initial parameters


func=bottleneck # set the function

############### Carry out optimization (same for any model) ########################
# Make extrapolation function:
func_ex = dadi.Numerics.make_extrap_log_func(func)
# perturb parameters
p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound,                                  lower_bound=lower_bound)
# optimize:
print('Beginning optimization ************************************************')
popt = dadi.Inference.optimize_log(p0, fs, func_ex, pts_l,
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=len(p0), maxiter=maxiter)
print('Finshed optimization **************************************************')

# Calculate the best-fit model AFS.
model = func_ex(popt, ns, pts_l)
# Likelihood of the data given the model AFS.
ll_model = dadi.Inference.ll_multinom(model, fs)
# calculate best fit theta
theta = dadi.Inference.optimal_sfs_scaling(model, fs)

###### model specific scaling of parameters (will depend on mu and L that you supply) #######
# param_names=("nuB","nuF","TF")
Nanc=theta / (4*mu*L)
nuB_scaled_dip=popt[0]*Nanc
nuF_scaled_dip=popt[1]*Nanc
TF_scaled_gen=popt[2]*2*Nanc
scaled_param_names=("Nanc_FromTheta_scaled_dip","nuB_scaled_dip","nuF_scaled_dip","T_scaled_gen")
scaled_popt=(Nanc,nuB_scaled_dip,nuF_scaled_dip,T_scaled_gen)


############### Write out output (same for any model) ########################
print('Writing out parameters **************************************************')

outputFile=open(str(outdir)+"/"+str(pop)+".dadi.inference."+str(modelName)+".runNum."+str(runNum)+"."+str(todaysdate)+".output","w")
# get all param names:
param_names_str='\t'.join(str(x) for x in param_names)
scaled_param_names_str='\t'.join(str(x) for x in scaled_param_names)
header=param_names_str+"\t"+scaled_param_names_str+"\ttheta\tLL\tmodelFunction\tmu\tL\tmaxiter\trunNumber\trundate\tinitialParameters\tupper_bound\tlower_bound" # add additional parameters theta, log-likelihood, model name, run number and rundate
popt_str='\t'.join(str(x) for x in popt) # get opt'd parameters as a tab-delim string
scaled_popt_str='\t'.join(str(x) for x in scaled_popt)
# joint together all the output fields, tab-separated:
output=[popt_str,scaled_popt_str,theta,ll_model,func.func_name,mu,L,maxiter,runNum,todaysdate,p0,upper_bound,lower_bound] # put all the output terms together
output='\t'.join(str(x) for x in output) # write out all the output fields
# this should result in a 2 row table that could be input into R / concatenated with other runs
outputFile.write(('{0}\n{1}\n').format(header,output))
outputFile.close()

############### Output SFS ########################
print('Writing out SFS **************************************************')

outputSFS=str(outdir)+"/"+str(pop)+".dadi.inference."+str(modelName)+".runNum."+str(runNum)+"."+str(todaysdate)+".expSFS"

# 20190117 -- fixed this to output EXPECTED sfs not obs sfs
model.to_file(outputSFS)
#outputSFS.close()

############### Output plot ########################
print('Making plots **************************************************')

#import pylab
import matplotlib.pyplot as plt
fig=plt.figure(1)
#pylab.ion()
outputFigure=str(str(outdir)+"/"+str(pop)+".dadi.inference."+str(modelName)+".runNum."+str(runNum)+"."+str(todaysdate)+".figure.png")
dadi.Plotting.plot_1d_comp_multinom(model, fs)
#pylab.show()
plt.savefig(outputFigure)
#pylab.clf()


###### exit #######
sys.exit()
