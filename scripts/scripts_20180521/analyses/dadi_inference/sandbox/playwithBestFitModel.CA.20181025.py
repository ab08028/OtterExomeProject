# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 16:25:06 2018

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


########## empirical sfs ##############
genotypeDate=20180806
sfsDate=20181019
sfs_dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/"+str(genotypeDate)+"/neutralSFS/"
pop="CA"
sfs=sfs_dir+"/"+pop+".unfolded.sfs.dadi.format."+str(sfsDate)+".txt"

fs=dadi.Spectrum.from_file(sfs) # this is unfolded
# fold the fs:
fs=fs.fold() # folded
ns = fs.sample_sizes # get sample size from SFS (in haploids)
pts_l = [ns[0]+5,ns[0]+15,ns[0]+25] # this should be slightly larger (+5) than sample size and increase by 10
########## Play with models #########

# figure out a way to automatically get best-fit params (from R?)

# try playing with the models

##################### Functions #########################

def one_bottleneck(params, ns, pts): 
    nuB,nuF,TB,TF = params
    xx = Numerics.default_grid(pts) # sets up grid
    phi = PhiManip.phi_1D(xx) # sets up initial phi for population 
    phi = Integration.one_pop(phi, xx, TB, nuB)  # bottleneck
    phi = Integration.one_pop(phi, xx, TF, nuF) # recovery 
    fs = Spectrum.from_phi(phi, ns, (xx,)) 
    return fs

def double_bottleneck(params, ns, pts): 
    nuB1,nuF1,nuB2,nuF2,TB1,TF1,TB2,TF2 = params
    xx = Numerics.default_grid(pts) # sets up grid
    phi = PhiManip.phi_1D(xx) # sets up initial phi for population 
    phi = Integration.one_pop(phi, xx, TB1, nuB1)  # bottleneck 1 
    phi = Integration.one_pop(phi, xx, TF1, nuF1) # recovery 1
    phi = Integration.one_pop(phi, xx, TB2, nuB2)  # bottleneck 2 
    phi = Integration.one_pop(phi, xx, TF2, nuF2) # recovery 2
    fs = Spectrum.from_phi(phi, ns, (xx,)) 
    return fs
    
func_one_bottleneck=one_bottleneck # set the function
func_double_bottleneck=double_bottleneck
############### Make extrapoloation functions ########################
# Make extrapolation function:
func_ex_one_bottleneck = dadi.Numerics.make_extrap_log_func(func_one_bottleneck)

func_ex_two_bottleneck = dadi.Numerics.make_extrap_log_func(func_double_bottleneck)

######################### explore one bottleneck model (best for CA) #######################
# from  /Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/20180806/CA/CA.bestModelRunParams.AIC.20181024.txt 
# eventually pull this out automatically
# original nuB,nuF,TB,TF=(0.00123692246704491, 94.2547453858811 ,0.000495116543844639, 0.00862901807609215)
# modify:
nuB,nuF,TB,TF=(0.01, 0.3 ,0.002, 0.002)
# get expected SFS under these parameters:
model = func_ex_one_bottleneck([nuB,nuF,TB,TF],ns,pts_l)

# calculate best fit theta
theta = dadi.Inference.optimal_sfs_scaling(model, fs)

# Likelihood of the data given the model AFS.
ll_model = dadi.Inference.ll_multinom(model, fs)


# plot fit:

######################### explore two bottleneck model (for CA) #######################
def double_bottleneck_fixed(params, ns, pts): 
    nuB1,nuF1,TB1,TF1 = params
    xx = Numerics.default_grid(pts) # sets up grid
    phi = PhiManip.phi_1D(xx) # sets up initial phi for population 
    phi = Integration.one_pop(phi, xx, TB1, nuB1)  # bottleneck 1 
    phi = Integration.one_pop(phi, xx, TF1, nuF1) # recovery 1
    phi = Integration.one_pop(phi, xx, 0.001, 0.01)  # bottleneck 2 
    phi = Integration.one_pop(phi, xx, 0.001, 0.1) # recovery 2
    fs = Spectrum.from_phi(phi, ns, (xx,)) 
    return fs
    
    
upper_bound = [10, 10, 10, 10]
lower_bound = [1e-4, 1e-4, 0, 0]
p0 = [0.01, 1, 0.0005, 0.01] # initial parameters
