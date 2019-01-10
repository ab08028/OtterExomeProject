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

# this model is the same as 1D.1Bottleneck but the duration of the bottleneck is fixed
modelName="1D.1Bottleneck.fixedDuration"

# for testing purposes:
sfs="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/20180806/neutralSFS/CA.unfolded.sfs.dadi.format.20181019.txt"
#sfs="/Users/annabelbeichman/Documents/UCLA/dadi/dadi/examples/YRI_CEU/YRI_CEU.fs" #testing
# project down

pop="CA"
mu=float(8.64411385098638e-09)
L=float(4193967)
############### Parse input arguments ########################
parser = argparse.ArgumentParser(description='Infer a '+ modelName +' model from a 1D folded SFS in dadi')
parser.add_argument("--runNum",required=True,help="iteration number (e.g. 1-50)")
parser.add_argument("--pop",required=True,help="population identifier, e.g. 'CA'")
parser.add_argument("--mu",required=True,help="supply mutation rate in mutation/bp/gen")
parser.add_argument("--L",required=True,help="number of called neutral sites that went into making SFS (monomorphic+polymorphic)")
parser.add_argument("--sfs",required=True,help="path to UNFOLDED SFS in dadi format (mask optional)")
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
fs=dadi.Spectrum.from_file(sfs) # this is unfolded
# fold the fs:
fs=fs.fold() # folded
############### Set up General Dadi Parameters ########################
ns = fs.sample_sizes # get sample size from SFS (in haploids)
pts_l = [ns[0]+5,ns[0]+15,ns[0]+25] # this should be slightly larger (+5) than sample size and increase by 10
############### Set up Specific Model -- this will change from script to script ########################

def bottleneck_fixedDur(params, ns, pts): 
    TB, TF = params # removed TB because that is fixed at fixed DUR
    xx = Numerics.default_grid(pts) # sets up grid
    phi = PhiManip.phi_1D(xx) # sets up initial phi for population 
    phi = Integration.one_pop(phi, xx, TB, 0.01)  # replace TB with fixed duration fors bottleneck
    phi = Integration.one_pop(phi, xx, TF, 0.1) # recovery 
    fs = Spectrum.from_phi(phi, ns, (xx,)) 
    return fs
param_names=("TB","TF")

def bottleneck(params, ns, pts): 
    nuB,nuF,TB,TF = params
    xx = Numerics.default_grid(pts) # sets up grid
    phi = PhiManip.phi_1D(xx) # sets up initial phi for population 
    phi = Integration.one_pop(phi, xx, TB, nuB)  # bottleneck
    phi = Integration.one_pop(phi, xx, TF, nuF) # recovery 
    fs = Spectrum.from_phi(phi, ns, (xx,)) 
    return fs
    
func=Demographics1D.bottlegrowth
param_names=("nuB","nuF","T")

upper_bound = [200, 200, 10 ]
lower_bound = [1e-4, 1e-4, 0]
p0 = [0.01,0.1,0.001] # initial parameters



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

############ for testing:
Nanc=theta/(4*mu*L)
nuB=popt[0]*Nanc
nuF=popt[1]*Nanc
TB=popt[2]*Nanc
TF=popt[3]*2*Nanc
Nanc,nuB,nuF,TB,TF,ll_model
############### Write out output (same for any model) ########################
print('Writing out parameters **************************************************')                                   

outputFile=open(str(outdir)+"/"+str(pop)+".dadi.inference."+str(modelName)+".runNum."+str(runNum)+"."+str(todaysdate)+".output","w")
# get all param names:
param_names_str='\t'.join(str(x) for x in param_names)
param_names_str=param_names_str+"\ttheta\tLL\tmodelFunction\tmu\tL\tmaxiter\trunNumber\trundate\tinitialParameters\tupper_bound\tlower_bound" # add additional parameters theta, log-likelihood, model name, run number and rundate
popt_str='\t'.join(str(x) for x in popt) # get opt'd parameters as a tab-delim string
# joint together all the output fields, tab-separated:
output=[popt_str,theta,ll_model,func.func_name,mu,L,maxiter,runNum,todaysdate,p0,upper_bound,lower_bound] # put all the output terms together
output='\t'.join(str(x) for x in output) # write out all the output fields
# this should result in a 2 row table that could be input into R / concatenated with other runs
outputFile.write(('{0}\n{1}\n').format(param_names_str,output))
outputFile.close()

############### Output SFS ########################
print('Writing out SFS **************************************************')                                   

outputSFS=str(outdir)+"/"+str(pop)+".dadi.inference."+str(modelName)+".runNum."+str(runNum)+"."+str(todaysdate)+".expSFS"

fs.to_file(outputSFS)
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
