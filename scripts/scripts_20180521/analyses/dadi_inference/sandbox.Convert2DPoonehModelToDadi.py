# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 16:13:58 2020

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
######### Want to compare Pooneh's model SFS to my empirical 2D sfs:
sfs="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/20181119/easySFS_projection/neutral/projection-20181221-hetFilter-0.75/dadi-plusMonomorphic/CA-AK.plusMonomorphic.sfs"
fs=dadi.Spectrum.from_file(sfs) # this is folded from easy SFS
# fold the fs:
#fs=fs.fold() # folded
# check if it's folded, if not folded, fold it
if fs.folded==False:
    fs=fs.fold()
else:
    fs=fs
ns = fs.sample_sizes # get sample size from SFS (in haploids)
pts_l = [ns[0]+5,ns[0]+15,ns[0]+25]
# want to compare Pooneh's model to this one
NancPooneh=float(3000)
# TDiv,nu1,nu2,m12,m21
paramsPooneh=[8000/(2*NancPooneh),1398/NancPooneh,3354/NancPooneh,4e-4*2*NancPooneh,4e-4*2*NancPooneh]
NancAnnabel=float(3800)
paramsAnnabel=[3927/(2*NancAnnabel),1828/NancAnnabel,2668/NancAnnabel,0.000158492*2*NancAnnabel,0.000158492*2*NancAnnabel]
def pooneh2DModel_ABMod((TDiv,nu1,nu2,m12,m21),ns,pts):
    xx = Numerics.default_grid(pts) # sets up grid
    phi = PhiManip.phi_1D(xx) # sets up initial phi for population 
    phi = PhiManip.phi_1D_to_2D(xx, phi)  # split into two pops
    phi = Integration.two_pops(phi, xx, TDiv, nu1, nu2, m12=m12, m21=m21)  # two pops at diff sizes with symmetric migration
    # allow for another size change; fix time at 35 gen
    # fix these at my bneck sizes (not pooneh's)
    phi = Integration.two_pops(phi, xx, 0.004, 0.05, 0.05, m12=0, m21=0) # fixing time at 35/(2*4000) to represent 35 gen ago and am fixing migraiton to be 0 during the bottleneck
    fs = Spectrum.from_phi(phi, ns, (xx,xx)) 
    return fs


func=pooneh2DModel_ABMod
modelName="pooneh2DModel_ABMod"
# wrap your function in Numerics.make_extrap_log_func to make extrapolation version of your function
func_ex = Numerics.make_extrap_log_func(func)

modelPooneh=func_ex(paramsPooneh,ns,pts_l) # this is relative to theta =1
ll_modelPooneh = dadi.Inference.ll_multinom(modelPooneh, fs)
ll_modelPooneh


modelAnnabel=func_ex(paramsAnnabel,ns,pts_l) # this is relative to theta =1
ll_modelAnnabel = dadi.Inference.ll_multinom(modelAnnabel, fs)
ll_modelAnnabel
# calculate best fit theta
thetaPooneh = dadi.Inference.optimal_sfs_scaling(modelPooneh, fs)
thetaAnnabel = dadi.Inference.optimal_sfs_scaling(modelAnnabel, fs)
thetaPooneh
thetaAnnabel

outdir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/20181119/comparingPoonehAnnabel2DModels/"
outputFigurePooneh=str(str(outdir)+"/"+"2D.expSFS.8KSplitTime.PoonehModel.figure.png")
outputFigureAnnabel=str(str(outdir)+"/"+"2D.expSFS.4KSplitTime.AnnabelModel.figure.png")
import matplotlib.pyplot as plt 

fig=plt.figure(1)
dadi.Plotting.plot_2d_comp_multinom(modelPooneh, fs)
pyplot.title((ll_modelPooneh))
plt.savefig(outputFigurePooneh)

fig=plt.figure(1)
dadi.Plotting.plot_2d_comp_multinom(modelAnnabel, fs)
pyplot.title((ll_modelAnnabel))
plt.savefig(outputFigureAnnabel)

