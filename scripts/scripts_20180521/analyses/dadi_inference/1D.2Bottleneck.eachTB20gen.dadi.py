# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 15:46:50 2018

@author: annabelbeichman
"""

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

modelName="1D.2Bottleneck"

############### Parse input arguments ########################
parser = argparse.ArgumentParser(description='Infer a '+ modelName +' model from a 1D folded SFS in dadi')
parser.add_argument("--runNum",required=True,help="iteration number (e.g. 1-50)")
parser.add_argument("--pop",required=True,help="population identifier, e.g. 'CA'")
parser.add_argument("--mu",required=True,help="supply mutation rate in mutation/bp/gen")
parser.add_argument("--L",required=True,help="number of called neutral sites that went into making SFS (monomorphic+polymorphic)")
parser.add_argument("--sfs",required=True,help="path to FOLDED SFS in dadi format from easySFS (mask optional)")
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
fs=dadi.Spectrum.from_file(sfs) # this is folded from easy sfs
# fold the fs:
#fs=fs.fold() # folded

############### Set up General Dadi Parameters ########################
ns = fs.sample_sizes # get sample size from SFS (in haploids)
pts_l = [ns[0]+5,ns[0]+15,ns[0]+25] # this should be slightly larger (+5) than sample size and increase by 10
############### Set up Specific Model -- this will change from script to script ########################
# Bernard says to fix bottleneck duration (TB1 and TB2) 
# 20 generations ~ 20/(3000*2); loosely setting it as 0.005 
def double_bottleneck(params, ns, pts): 
    nuB1,nuF1,nuB2,nuF2,TB1,TF1,TB2,TF2 = params
    xx = Numerics.default_grid(pts) # sets up grid
    phi = PhiManip.phi_1D(xx) # sets up initial phi for population 
    phi = Integration.one_pop(phi, xx, 0.005, nuB1)  # bottleneck 1 
    phi = Integration.one_pop(phi, xx, TF1, nuF1) # recovery 1
    phi = Integration.one_pop(phi, xx, 0.005, nuB2)  # bottleneck 2 
    phi = Integration.one_pop(phi, xx, TF2, nuF2) # recovery 2
    fs = Spectrum.from_phi(phi, ns, (xx,)) 
    return fs
param_names=("nuB1","nuF1","nuB2","nuF2","TF1","TF2")

# 20181024 changing upper bound on pop sizes to 10 because know it doesn't grow up to 100* Nanc; and lowering lower bounds to 1e-4 ; see what happens
# 20181031 changed bounds based on bernard suggestions; went with .1 as upper bound on T intervals (similar to humans)
# and 1e-5 as lower bound on times (0 lets intervals collapse which can screw up optimization)

upper_bound = [10, 10, 10, 10, .1, .1, .1, .1]
lower_bound = [1e-4, 1e-4, 1e-4, 1e-4, 1e-5, 1e-5, 1e-5, 1e-5]
p0 = [0.01, 1, 0.01, 0.3, 0.0005,0.01,0.0005,0.0005] # initial parameters


func=double_bottleneck # set the function

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
# param_names=("nuB1","nuF1","nuB2","nuF2","TF1","TF2")
Nanc=theta / (4*mu*L)
nuB1_scaled_dip=popt[0]*Nanc
nuF1_scaled_dip=popt[1]*Nanc
nuB2_scaled_dip=popt[2]*Nanc
nuF2_scaled_dip=popt[3]*Nanc
TF1_scaled_gen=popt[4]*2*Nanc
TF2_scaled_gen=popt[5]*2*Nanc
scaled_param_names=("Nanc_FromTheta_scaled_dip","nuB1_scaled_dip","nuF1_scaled_dip","nuB2_scaled_dip","nuF2_scaled_dip","TF1_scaled_gen","TF2_scaled_gen")
scaled_popt=(Nanc,nuB1_scaled_dip,nuF1_scaled_dip,nuB2_scaled_dip,nuF2_scaled_dip,TF1_scaled_gen,TF2_scaled_gen)


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
