# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 15:32:23 2018

@author: annabelbeichman
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 16:46:27 2018

@author: annabelbeichman
"""
import matplotlib
matplotlib.use('Agg') # so it doesn't pop up graphics on hoffman
import sys
import argparse
import dadi
from dadi import Numerics, PhiManip, Integration, Spectrum, Demographics1D
# note: module  Demographics1D has the following models file:///Users/annabelbeichman/Documents/UCLA/dadi/dadi/doc/api/dadi.Demographics1D-module.html

from numpy import array # don't comment this out
import datetime
todaysdate=datetime.datetime.today().strftime('%Y%m%d')

modelName="1D.1Epoch"
############### Parse input arguments ########################
parser = argparse.ArgumentParser(description='Infer a '+ modelName +' model from a 1D folded SFS in dadi')
parser.add_argument("--runNum",required=True,help="iteration number (e.g. 1-50)")
parser.add_argument("--pop",required=True,help="population identifier, e.g. 'CA'")
parser.add_argument("--mu",required=True,help="supply mutation rate in mutation/bp/gen")
parser.add_argument("--L",required=True,help="number of called neutral sites that went into making SFS (monomorphic+polymorphic)")
parser.add_argument("--sfs",required=True,help="path to FOLDED SFS in dadi format from easysfs (mask optional)")
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
#testfile="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/20181119/easySFS_projection/neutral/projection-20181221-hetFilter-0.75/dadi-plusMonomorphic/AK-15.plusMonomorphic.sfs"
fs=dadi.Spectrum.from_file(sfs) # this is folded if from easy SFS
# fold the fs:
#fs=fs.fold() # folded

############### Set up General Dadi Parameters ########################
ns = fs.sample_sizes # get sample size from SFS (in haploids)
pts_l = [ns[0]+5,ns[0]+15,ns[0]+25] # this should be slightly larger (+5) than sample size and increase by 10
############### Set up Specific Model -- this will change from script to script ########################
func = Demographics1D.snm # standard neutral model (?) #
# about this model "The model with no size changes is dadi.Demographics1D.snm. It has no free parameters, so no optimization to be done, just calculating theta.
#(Backstory: The integration always starts with the frequency distribution for a population that has stable size for a very long time. Any integration you do beyond that is introducing demographic events.)" https://groups.google.com/forum/#!msg/dadi-user/d0CRBUOO874/frLIBLE2FQAJ

# don't want to do parameter optimization

############### Carry out optimization (same for any model) ########################
# Make extrapolation function: (still need to do this for SNM)
# it just won't take any parameters
func_ex = dadi.Numerics.make_extrap_log_func(func)
# there are no parameters; don't need to optimize
                          
# Calculate the best-fit model AFS.
# first argument is "unused" -- can be anything, but isn't popt
# then you still use ns and pts_l
# https://groups.google.com/forum/#!searchin/dadi-user/snm%7Csort:date/dadi-user/ZEpNEt6zU9s/tcCvyiIjBQAJ
model = func_ex("N/A", ns, pts_l)
# Likelihood of the data given the model AFS.
ll_model = dadi.Inference.ll_multinom(model, fs)
# calculate best fit theta
theta = dadi.Inference.optimal_sfs_scaling(model, fs)

###### model specific scaling of parameters (will depend on mu and L that you supply) #######
#param_names= ("nu","T")

Nanc=theta / (4*mu*L)

scaled_param_names=("Nanc_FromTheta_scaled_dip")
scaled_popt=(Nanc)


############### Write out output (same for any model) ########################
print('Writing out parameters **************************************************')                                   

outputFile=open(str(outdir)+"/"+str(pop)+".dadi.inference."+str(modelName)+".runNum."+str(runNum)+"."+str(todaysdate)+".output","w")
# get all param names:
#param_names_str='\t'.join(str(x) for x in param_names)
scaled_param_names_str='\t'.join(str(x) for x in scaled_param_names)
header=scaled_param_names_str+"\ttheta\tLL\tmodelFunction\tmu\tL\tmaxiter\trunNumber\trundate\tinitialParameters\tupper_bound\tlower_bound" # add additional parameters theta, log-likelihood, model name, run number and rundate
scaled_popt_str='\t'.join(str(x) for x in scaled_popt)
# joint together all the output fields, tab-separated:
output=[scaled_popt_str,theta,ll_model,func.func_name,mu,L,maxiter,runNum,todaysdate] # put all the output terms together
output='\t'.join(str(x) for x in output) # write out all the output fields
# this should result in a 2 row table that could be input into R / concatenated with other runs
outputFile.write(('{0}\n{1}\n').format(header,output))
outputFile.close()


############### Output SFS ########################
print('Writing out SFS **************************************************')                                   

outputSFS=str(outdir)+"/"+str(pop)+".dadi.inference."+str(modelName)+".runNum."+str(runNum)+"."+str(todaysdate)+".expSFS"

# 20190117 -- fixed this! used to write out fs (obs sfs) but needed to write out exp sfs (caught it in time)
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
