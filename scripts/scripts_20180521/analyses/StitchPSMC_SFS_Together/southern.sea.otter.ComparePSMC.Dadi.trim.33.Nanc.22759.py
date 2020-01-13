# dadi function for sso msmc model; trim point = 33; Nanc=22759
# note: order of variables goes from oldest --> youngest; so nu0 is Nanc, T0 is time for Nanc. nu32 would be most recent if there were 33 time indices.
# This is reverse order for time index numbering in MSMC so don't get confused. They are listed in descending order because that's how they come from MSMC.
import dadi
from dadi import Numerics, PhiManip, Integration, Spectrum, Plotting, Inference
import numpy as np
import matplotlib.pyplot as pyplot

################# set up dirs #####################
pop="CA"
sfsdir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/20181119/easySFS_projection/neutral/projection-20181221-hetFilter-0.75/dadi-plusMonomorphic/"
outdir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/StitchPSMC_SFS_Together/"+pop # outdir
modelName="sso_model_trim_33"
############### Input data ####################################
sfs=sfsdir+"CA-13.plusMonomorphic.sfs" # get CA SFS 
fs=dadi.Spectrum.from_file(sfs) # this is folded if from easy SFS
# fold the fs:
#fs=fs.fold() # folded
# check if it's folded, if not folded, fold it
if fs.folded==False:
    fs=fs.fold()
else:
    fs=fs
############### Set up General Dadi Parameters ########################
ns = fs.sample_sizes # get sample size from SFS (in haploids)
pts_l = [ns[0]+5,ns[0]+15,ns[0]+25] # this should be slightly larger (+5) than sample size and increase by 10

Nanc = 22759
# params in dadi units (nu = Nx/Nanc; T = gen/(2*Nanc)) in same order that they should go into function (all nus from nu32 --> nu0 followed by Ts from T32 --> T0) -- note that order is from recent --> old because that is how it comes out of msmc. But then dadi function starts with oldest (nu0) first.
params = (0.0881030544009925, 0.0339730033926217, 0.0929355795355233, 0.146484994149822, 0.185969764821238, 0.218631178707224, 0.244857652102702, 0.262357053342665, 0.270286941253667, 0.269710116013023, 0.257097104104347, 0.257097104104347, 0.232676303911964, 0.232676303911964, 0.208900140554491, 0.208900140554491, 0.191032772098617, 0.191032772098617, 0.18040432146964, 0.18040432146964, 0.177336636081359, 0.177336636081359, 0.182994441404337, 0.182994441404337, 0.200959926621754, 0.200959926621754, 0.240224202954715, 0.240224202954715, 0.323417652253806, 0.323417652253806, 0.512531459731541, 0.512531459731541, 1, 0.00723549801, 0.00742344193499999, 0.00762142478500003, 0.00783031066999984, 0.00805070955000007, 0.00828427340000002, 0.00853156135000007, 0.00879409829999996, 0.00907340914999987, 0.00937101879999997, 0.00968845215000001, 0.0100287589999999, 0.0103932101, 0.0107858718499999, 0.0112085233, 0.0116665016000001, 0.0121628565499998, 0.0127036877500001, 0.01329484065, 0.0139436855999999, 0.0146588637000001, 0.0154518116999998, 0.0163349829500007, 0.0173256596499985, 0.0184436655000015, 0.0197169569999992, 0.0211808609999998, 0.0228735000000012, 0.0248685774999983, 0.0272372555000004, 0.0301116919999995, 0.0336596260000006, 0.0381606225000007)


################### copy/paste MSMC information from 
def sso_model_trim_33((nu32, nu31, nu30, nu29, nu28, nu27, nu26, nu25, nu24, nu23, nu22, nu21, nu20, nu19, nu18, nu17, nu16, nu15, nu14, nu13, nu12, nu11, nu10, nu9, nu8, nu7, nu6, nu5, nu4, nu3, nu2, nu1, nu0, T32, T31, T30, T29, T28, T27, T26, T25, T24, T23, T22, T21, T20, T19, T18, T17, T16, T15, T14, T13, T12, T11, T10, T9, T8, T7, T6, T5, T4, T3, T2, T1, T0),ns,pts):
	xx = Numerics.default_grid(pts)
	# intialize phi with ancestral pop: (nu0 = 1)
	phi = PhiManip.phi_1D(xx,nu=nu0)
	# stays at nu0 for T0 duration of time:
	# followed by a number of time steps, with associated pop changes:
	phi = Integration.one_pop(phi, xx, T0, nu0)
	phi = Integration.one_pop(phi, xx, T1, nu1)
	phi = Integration.one_pop(phi, xx, T2, nu2)
	phi = Integration.one_pop(phi, xx, T3, nu3)
	phi = Integration.one_pop(phi, xx, T4, nu4)
	phi = Integration.one_pop(phi, xx, T5, nu5)
	phi = Integration.one_pop(phi, xx, T6, nu6)
	phi = Integration.one_pop(phi, xx, T7, nu7)
	phi = Integration.one_pop(phi, xx, T8, nu8)
	phi = Integration.one_pop(phi, xx, T9, nu9)
	phi = Integration.one_pop(phi, xx, T10, nu10)
	phi = Integration.one_pop(phi, xx, T11, nu11)
	phi = Integration.one_pop(phi, xx, T12, nu12)
	phi = Integration.one_pop(phi, xx, T13, nu13)
	phi = Integration.one_pop(phi, xx, T14, nu14)
	phi = Integration.one_pop(phi, xx, T15, nu15)
	phi = Integration.one_pop(phi, xx, T16, nu16)
	phi = Integration.one_pop(phi, xx, T17, nu17)
	phi = Integration.one_pop(phi, xx, T18, nu18)
	phi = Integration.one_pop(phi, xx, T19, nu19)
	phi = Integration.one_pop(phi, xx, T20, nu20)
	phi = Integration.one_pop(phi, xx, T21, nu21)
	phi = Integration.one_pop(phi, xx, T22, nu22)
	phi = Integration.one_pop(phi, xx, T23, nu23)
	phi = Integration.one_pop(phi, xx, T24, nu24)
	phi = Integration.one_pop(phi, xx, T25, nu25)
	phi = Integration.one_pop(phi, xx, T26, nu26)
	phi = Integration.one_pop(phi, xx, T27, nu27)
	phi = Integration.one_pop(phi, xx, T28, nu28)
	phi = Integration.one_pop(phi, xx, T29, nu29)
	phi = Integration.one_pop(phi, xx, T30, nu30)
	phi = Integration.one_pop(phi, xx, T31, nu31)
	phi = Integration.one_pop(phi, xx, T32, nu32)
	# get expected SFS:
	fs = Spectrum.from_phi(phi,ns,(xx,))
	return fs


# wrap your function in Numerics.make_extrap_log_func to make extrapolation version of your function
extrap_sso_model_trim_33_function = Numerics.make_extrap_log_func(sso_model_trim_33)


msmc_model_expSFS=extrap_sso_model_trim_33_function(params,ns,pts_l)


ll_msmc_model = dadi.Inference.ll_multinom(msmc_model_expSFS, fs)

# plot an image:
#import pylab
import matplotlib.pyplot as plt 
#fig=plt.figure(1)
#pylab.ion()
outputFigure=str(str(outdir)+"/"+str(pop)+"."+str(modelName)+".wDip.figure.png")
dadi.Plotting.plot_1d_comp_multinom(msmc_model_expSFS, fs)
pyplot.title(str(modelName)+" LL = "+str(ll_msmc_model))
plt.savefig(outputFigure)

###################### Try without that semi-recent dip ###############
# converted that nu31 interval that was 0.0339730033926217 to  0.0929355795355233 (same as previous interval)
params_noDip = (0.0881030544009925, 0.0929355795355233, 0.0929355795355233, 0.146484994149822, 0.185969764821238, 0.218631178707224, 0.244857652102702, 0.262357053342665, 0.270286941253667, 0.269710116013023, 0.257097104104347, 0.257097104104347, 0.232676303911964, 0.232676303911964, 0.208900140554491, 0.208900140554491, 0.191032772098617, 0.191032772098617, 0.18040432146964, 0.18040432146964, 0.177336636081359, 0.177336636081359, 0.182994441404337, 0.182994441404337, 0.200959926621754, 0.200959926621754, 0.240224202954715, 0.240224202954715, 0.323417652253806, 0.323417652253806, 0.512531459731541, 0.512531459731541, 1, 0.00723549801, 0.00742344193499999, 0.00762142478500003, 0.00783031066999984, 0.00805070955000007, 0.00828427340000002, 0.00853156135000007, 0.00879409829999996, 0.00907340914999987, 0.00937101879999997, 0.00968845215000001, 0.0100287589999999, 0.0103932101, 0.0107858718499999, 0.0112085233, 0.0116665016000001, 0.0121628565499998, 0.0127036877500001, 0.01329484065, 0.0139436855999999, 0.0146588637000001, 0.0154518116999998, 0.0163349829500007, 0.0173256596499985, 0.0184436655000015, 0.0197169569999992, 0.0211808609999998, 0.0228735000000012, 0.0248685774999983, 0.0272372555000004, 0.0301116919999995, 0.0336596260000006, 0.0381606225000007)

msmc_model_expSFS_noDip=extrap_sso_model_trim_33_function(params_noDip,ns,pts_l)


ll_msmc_model_noDip = dadi.Inference.ll_multinom(msmc_model_expSFS_noDip, fs)

outputFigure=str(str(outdir)+"/"+str(pop)+"."+str(modelName)+".noDIP.figure.png")
dadi.Plotting.plot_1d_comp_multinom(msmc_model_expSFS_noDip, fs)
pyplot.title(str(modelName)+".noDip LL = "+str(ll_msmc_model_noDip))
plt.savefig(outputFigure)

########################## Try to add the contraction ###########
# how to appropriately scale?     

contractionGen=35 # 35 gen
contractionSize=195 # diploids 
# rescale by nanc:
contractionGen_RescBy2Nanc = float(contractionGen)/(2*Nanc)
contractionSize_RescByNanc = float(contractionSize)/Nanc
## adding them to end of parameter vector:
def sso_model_trim_33_plusContraction((nu32, nu31, nu30, nu29, nu28, nu27, nu26, nu25, nu24, nu23, nu22, nu21, nu20, nu19, nu18, nu17, nu16, nu15, nu14, nu13, nu12, nu11, nu10, nu9, nu8, nu7, nu6, nu5, nu4, nu3, nu2, nu1, nu0, T32, T31, T30, T29, T28, T27, T26, T25, T24, T23, T22, T21, T20, T19, T18, T17, T16, T15, T14, T13, T12, T11, T10, T9, T8, T7, T6, T5, T4, T3, T2, T1, T0,contractionSize_RescByNanc,contractionGen_RescBy2Nanc),ns,pts):
	xx = Numerics.default_grid(pts)
	# intialize phi with ancestral pop: (nu0 = 1)
	phi = PhiManip.phi_1D(xx,nu=nu0)
	# stays at nu0 for T0 duration of time:
	# followed by a number of time steps, with associated pop changes:
	phi = Integration.one_pop(phi, xx, T0, nu0)
	phi = Integration.one_pop(phi, xx, T1, nu1)
	phi = Integration.one_pop(phi, xx, T2, nu2)
	phi = Integration.one_pop(phi, xx, T3, nu3)
	phi = Integration.one_pop(phi, xx, T4, nu4)
	phi = Integration.one_pop(phi, xx, T5, nu5)
	phi = Integration.one_pop(phi, xx, T6, nu6)
	phi = Integration.one_pop(phi, xx, T7, nu7)
	phi = Integration.one_pop(phi, xx, T8, nu8)
	phi = Integration.one_pop(phi, xx, T9, nu9)
	phi = Integration.one_pop(phi, xx, T10, nu10)
	phi = Integration.one_pop(phi, xx, T11, nu11)
	phi = Integration.one_pop(phi, xx, T12, nu12)
	phi = Integration.one_pop(phi, xx, T13, nu13)
	phi = Integration.one_pop(phi, xx, T14, nu14)
	phi = Integration.one_pop(phi, xx, T15, nu15)
	phi = Integration.one_pop(phi, xx, T16, nu16)
	phi = Integration.one_pop(phi, xx, T17, nu17)
	phi = Integration.one_pop(phi, xx, T18, nu18)
	phi = Integration.one_pop(phi, xx, T19, nu19)
	phi = Integration.one_pop(phi, xx, T20, nu20)
	phi = Integration.one_pop(phi, xx, T21, nu21)
	phi = Integration.one_pop(phi, xx, T22, nu22)
	phi = Integration.one_pop(phi, xx, T23, nu23)
	phi = Integration.one_pop(phi, xx, T24, nu24)
	phi = Integration.one_pop(phi, xx, T25, nu25)
	phi = Integration.one_pop(phi, xx, T26, nu26)
	phi = Integration.one_pop(phi, xx, T27, nu27)
	phi = Integration.one_pop(phi, xx, T28, nu28)
	phi = Integration.one_pop(phi, xx, T29, nu29)
	phi = Integration.one_pop(phi, xx, T30, nu30)
	phi = Integration.one_pop(phi, xx, T31, nu31)
	phi = Integration.one_pop(phi, xx, T32, nu32)
	### add contraction:
	phi = Integration.one_pop(phi, xx,contractionGen_RescBy2Nanc,  contractionSize_RescByNanc)
	# get expected SFS:
	fs = Spectrum.from_phi(phi,ns,(xx,))
	return fs



# wrap your function in Numerics.make_extrap_log_func to make extrapolation version of your function
extrap_sso_model_trim_33_plusContraction_function = Numerics.make_extrap_log_func(sso_model_trim_33_plusContraction)

params_plusContraction=(0.0881030544009925, 0.0339730033926217, 0.0929355795355233, 0.146484994149822, 0.185969764821238, 0.218631178707224, 0.244857652102702, 0.262357053342665, 0.270286941253667, 0.269710116013023, 0.257097104104347, 0.257097104104347, 0.232676303911964, 0.232676303911964, 0.208900140554491, 0.208900140554491, 0.191032772098617, 0.191032772098617, 0.18040432146964, 0.18040432146964, 0.177336636081359, 0.177336636081359, 0.182994441404337, 0.182994441404337, 0.200959926621754, 0.200959926621754, 0.240224202954715, 0.240224202954715, 0.323417652253806, 0.323417652253806, 0.512531459731541, 0.512531459731541, 1, 0.00723549801, 0.00742344193499999, 0.00762142478500003, 0.00783031066999984, 0.00805070955000007, 0.00828427340000002, 0.00853156135000007, 0.00879409829999996, 0.00907340914999987, 0.00937101879999997, 0.00968845215000001, 0.0100287589999999, 0.0103932101, 0.0107858718499999, 0.0112085233, 0.0116665016000001, 0.0121628565499998, 0.0127036877500001, 0.01329484065, 0.0139436855999999, 0.0146588637000001, 0.0154518116999998, 0.0163349829500007, 0.0173256596499985, 0.0184436655000015, 0.0197169569999992, 0.0211808609999998, 0.0228735000000012, 0.0248685774999983, 0.0272372555000004, 0.0301116919999995, 0.0336596260000006, 0.0381606225000007,contractionSize_RescByNanc,contractionGen_RescBy2Nanc)

msmc_model_expSFS_wDip_wContraction=extrap_sso_model_trim_33_plusContraction_function(params_plusContraction,ns,pts_l)


ll_msmc_model_wDip_wContraction = dadi.Inference.ll_multinom(msmc_model_expSFS_wDip_wContraction, fs)

ll_msmc_model_wDip_wContraction

outputFigure=str(str(outdir)+"/"+str(pop)+"."+str(modelName)+".wDIP.AddRecentContraction.figure.png")
dadi.Plotting.plot_1d_comp_multinom(msmc_model_expSFS_wDip_wContraction, fs)
pyplot.title(str(modelName)+".wRecentContraction LL = "+str(ll_msmc_model_wDip_wContraction))
plt.savefig(outputFigure)

########################## try with No dip and dadi contraction ######
# converted that nu31 interval that was 0.0339730033926217 to  0.0929355795355233 (same as previous interval)
params_plusContraction_noDIP= (0.0881030544009925, 0.0929355795355233, 0.0929355795355233, 0.146484994149822, 0.185969764821238, 0.218631178707224, 0.244857652102702, 0.262357053342665, 0.270286941253667, 0.269710116013023, 0.257097104104347, 0.257097104104347, 0.232676303911964, 0.232676303911964, 0.208900140554491, 0.208900140554491, 0.191032772098617, 0.191032772098617, 0.18040432146964, 0.18040432146964, 0.177336636081359, 0.177336636081359, 0.182994441404337, 0.182994441404337, 0.200959926621754, 0.200959926621754, 0.240224202954715, 0.240224202954715, 0.323417652253806, 0.323417652253806, 0.512531459731541, 0.512531459731541, 1, 0.00723549801, 0.00742344193499999, 0.00762142478500003, 0.00783031066999984, 0.00805070955000007, 0.00828427340000002, 0.00853156135000007, 0.00879409829999996, 0.00907340914999987, 0.00937101879999997, 0.00968845215000001, 0.0100287589999999, 0.0103932101, 0.0107858718499999, 0.0112085233, 0.0116665016000001, 0.0121628565499998, 0.0127036877500001, 0.01329484065, 0.0139436855999999, 0.0146588637000001, 0.0154518116999998, 0.0163349829500007, 0.0173256596499985, 0.0184436655000015, 0.0197169569999992, 0.0211808609999998, 0.0228735000000012, 0.0248685774999983, 0.0272372555000004, 0.0301116919999995, 0.0336596260000006, 0.0381606225000007,contractionSize_RescByNanc,contractionGen_RescBy2Nanc)

msmc_model_expSFS_wContraction_noDIP=extrap_sso_model_trim_33_plusContraction_function(params_plusContraction_noDIP,ns,pts_l)


ll_msmc_model_wContraction_noDIP = dadi.Inference.ll_multinom(msmc_model_expSFS_wContraction_noDIP, fs)

ll_msmc_model_wContraction_noDIP

outputFigure=str(str(outdir)+"/"+str(pop)+"."+str(modelName)+".noDIP.AddRecentContraction.figure.png")
dadi.Plotting.plot_1d_comp_multinom(msmc_model_expSFS_wContraction_noDIP, fs)
pyplot.title(str(modelName)+".noDIP.wRecentContraction LL = "+str(ll_msmc_model_wContraction_noDIP))
plt.savefig(outputFigure)