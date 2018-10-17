# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 10:40:25 2018

@author: annabelbeichman
"""
import sys
import dadi
from dadi import Numerics, PhiManip, Integration, Spectrum
from numpy import array
import datetime
todaysdate=datetime.datetime.today().strftime('%Y%m%d')

# Designing my 1D otter inference models

pop="CA"
fs=dadi.Spectrum.from_file("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/20180806/neutralSFS/CA_all_9_rmAllHet_rmRelativesAdmixed_passingAllFilters_allCalled.filtered.sfs.20181009.dadi.format.out")
ns = fs.sample_sizes # 34 (should this be 18 for folded? -- ask Bernard)
L=4193488 # for California: total called neutral sites # get this from data table at some point
mu=8.64411385098638e-09 # mu I got in my otter paper
# Implementing single bottleneck model:
# try this tutorial: https://bitbucket.org/gutenkunstlab/dadi/src/476f25172d66fcbdf65940d7bbe3a6c0e6d6013c/examples/YRI_CEU/YRI_CEU.py?at=master&fileviewer=file-view-default


pts_l = [40,50,60] # arbitrary
# good results are obtained by setting the smallest grid size slightly larger than the largest population sample size (haploids). So my haploid size is 34, so 40,50,60 is okay for now (check this)
# from dadi manual:
def bottleneck(params, ns, pts): 
    nuB,nuF,TB,TF = params
    xx = Numerics.default_grid(pts) # sets up grid
    phi = PhiManip.phi_1D(xx) # sets up initial phi for population 
    phi = Integration.one_pop(phi, xx, TB, nuB)  # bottleneck
    phi = Integration.one_pop(phi, xx, TF, nuF) # recovery 
    fs = Spectrum.from_phi(phi, ns, (xx,)) 
    return fs
#param_names=("nuB: bottleneck size","nuF: recovery size","TB: bottleneck duration","TF: recovery duration")
param_names=("nuB","nuF","TB","TF")
# this model says that at time TF+TB before present, the population goes into a bottleneck of size nuB
# and then at time TF before present it recovers to size nuF

func=bottleneck

# Parameters are: (nuB (bneck size),nuF (recovery size),TB,TF)
# sizes are relative to Nref (ancestral size)
# In our analyses, we often set the upper bound on times to be 10 and the upper bound on migration rates to be 20.(dadi manual)
upper_bound = [100, 100, 10, 10]
# we often set the lower bound on population sizes to be 1e- 2 or 1e-3 
lower_bound = [1e-3, 1e-3, 0, 0]

# This is our initial guess for the parameters, which is somewhat arbitrary.
# picking 0.01 for bneck size (100/10,000) and 0.3 for recovery size (3000/10,000) and 1 and 1 or time intervals (arbitrary)

p0 = [0.01,0.3,0.001,0.001] # made initial times wayyy smaller (from 1 -> 0.005) and it helps a lot (much more recent/short term events)

# Make the extrapolating version of our demographic model function.
func_ex = dadi.Numerics.make_extrap_log_func(func)


# Perturb our parameters before optimization. This does so by taking each
# parameter a up to a factor of two up or down. randomly perturbs 
p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound,                                  lower_bound=lower_bound)
# is this enough for changing starting parameters?
# Do the optimization. By default we assume that theta is a free parameter,
# since it's trivial to find given the other parameters. If you want to fix
# theta, add a multinom=False to the call.
# The maxiter argument restricts how long the optimizer will run. 
# ****** For real 
# runs, you will want to set this value higher (at least 10), *****
#  to encourage
# better convergence. ***  You will also want to run optimization several times
# using multiple sets of intial parameters, ****  to be confident you've actually
# found the true maximum likelihood parameters.

print('Beginning optimization ************************************************')
popt = dadi.Inference.optimize_log(p0, fs, func_ex, pts_l, 
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=len(p0), maxiter=10)
                                   ### eventually change maxiter to 10+
# The verbose argument controls how often progress of the optimizer should be
# printed. It's useful to keep track of optimization process.
print('Finshed optimization **************************************************')                                   
# Calculate the best-fit model AFS.
model = func_ex(popt, ns, pts_l)
# Likelihood of the data given the model AFS.
ll_model = dadi.Inference.ll_multinom(model, fs)
# calculate best fit theta
theta = dadi.Inference.optimal_sfs_scaling(model, fs)
######## sandbox #########
popt_str='\t'.join(str(x) for x in popt)
output=[popt_str,theta,ll_model,func.func_name,todaysdate]
print('\t'.join(str(x) for x in param_names)),
print('theta\tLL\tmodel\trundate')
print('\n')
print('\t'.join(str(x) for x in output))


##################################
print('Model fit: {0}'.format(func.func_name))
print('Parameters fit: {0}'.format(param_names))
print('Best-fit parameters: {0}'.format(popt))
print('Best-fit parameters with names: {0}'.format(zip(param_names,popt)))
print('Note units: sizes (nu) are relative to Nanc (diploid); times are in 2*Nanc generations')
# Calculate the best-fit model AFS.
model = func_ex(popt, ns, pts_l)
# Likelihood of the data given the model AFS.
ll_model = dadi.Inference.ll_multinom(model, fs)
print('Maximum log composite likelihood: {0}'.format(ll_model))
theta = dadi.Inference.optimal_sfs_scaling(model, fs)
print('Optimal value of theta: {0}'.format(theta))
# calculate Nref
Nref=theta/(4*mu*L) # diploids!
print('Optimal value of Nref based on theta ({0}), mu ({1}) and L ({2}): Nref = {3} diploids '.format(theta, mu,L,Nref))
####### Note that Nref is in diploids! 
# Plot:
import pylab
#pylab.figure(1)
dadi.Plotting.plot_1d_comp_multinom(model, fs)
# save figure:
#pylab.show()
#pylab.show()
pylab.savefig('CA.model1.png', dpi=50)

## make this script nicer to that I write out the results

sys.exit()