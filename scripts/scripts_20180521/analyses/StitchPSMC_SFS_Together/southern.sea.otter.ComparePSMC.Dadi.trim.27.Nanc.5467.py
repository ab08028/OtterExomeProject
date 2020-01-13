### want to better automate these scripts for different trim values (can do in R easily enough I think)

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
modelName="sso_model_trim_27" # change by hand for now
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

################### Copy in parameters from the model file you output in R ############

from dadi import Numerics, PhiManip, Integration, Spectrum, Plotting, Inference
import numpy as np
import matplotlib.pyplot as pyplot


Nanc = 5467
# params in dadi units (nu = Nx/Nanc; T = gen/(2*Nanc)) in same order that they should go into function (all nus from nu32 --> nu0 followed by Ts from T32 --> T0) -- note that order is from recent --> old because that is how it comes out of msmc. But then dadi function starts with oldest (nu0) first.
params = (0.36675344664418, 0.141422067280315, 0.386870175412936, 0.609784494435127, 0.774150824662306, 0.91011303614748, 1.01928801965413, 1.09213414017289, 1.12514450221579, 1.12274330685933, 1.07023813979649, 1.07023813979649, 0.968579772770968, 0.968579772770968, 0.869604885707006, 0.869604885707006, 0.795226999398677, 0.795226999398677, 0.750983120146509, 0.750983120146509, 0.738213027247672, 0.738213027247672, 0.761765214135537, 0.761765214135537, 0.836551538729166, 0.836551538729166, 1, 0.030119771118, 0.030902139933, 0.0317262985630002, 0.0325958441059994, 0.0335133156900004, 0.0344855901200002, 0.0355149949300004, 0.0366078779399999, 0.0377705869699996, 0.03900946984, 0.0403308743700001, 0.0417474961999999, 0.0432646251800001, 0.0448991888299998, 0.04665859294, 0.0485650548800005, 0.0506312702899994, 0.0528826304500005, 0.05534346867, 0.0580444660799995, 0.0610215936600004, 0.0643224600599991, 0.0679989058100032, 0.072122872869994, 0.0767768829000064, 0.0820773125999968, 0.0881712197999992)



def sso_model_trim_27((nu26, nu25, nu24, nu23, nu22, nu21, nu20, nu19, nu18, nu17, nu16, nu15, nu14, nu13, nu12, nu11, nu10, nu9, nu8, nu7, nu6, nu5, nu4, nu3, nu2, nu1, nu0, T26, T25, T24, T23, T22, T21, T20, T19, T18, T17, T16, T15, T14, T13, T12, T11, T10, T9, T8, T7, T6, T5, T4, T3, T2, T1, T0),ns,pts):
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
	# get expected SFS:
	fs = Spectrum.from_phi(phi,ns,(xx,))
	return fs


# wrap your function in Numerics.make_extrap_log_func to make extrapolation version of your function
extrap_sso_model_trim_27_function = Numerics.make_extrap_log_func(sso_model_trim_27)


msmc_model_expSFS=extrap_sso_model_trim_27_function(params,ns,pts_l)


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
# converted that nu31 interval that was 0.141422067280315 to  0.386870175412936 (same as previous interval)
params_noDip = (0.36675344664418, 0.386870175412936, 0.386870175412936, 0.609784494435127, 0.774150824662306, 0.91011303614748, 1.01928801965413, 1.09213414017289, 1.12514450221579, 1.12274330685933, 1.07023813979649, 1.07023813979649, 0.968579772770968, 0.968579772770968, 0.869604885707006, 0.869604885707006, 0.795226999398677, 0.795226999398677, 0.750983120146509, 0.750983120146509, 0.738213027247672, 0.738213027247672, 0.761765214135537, 0.761765214135537, 0.836551538729166, 0.836551538729166, 1, 0.030119771118, 0.030902139933, 0.0317262985630002, 0.0325958441059994, 0.0335133156900004, 0.0344855901200002, 0.0355149949300004, 0.0366078779399999, 0.0377705869699996, 0.03900946984, 0.0403308743700001, 0.0417474961999999, 0.0432646251800001, 0.0448991888299998, 0.04665859294, 0.0485650548800005, 0.0506312702899994, 0.0528826304500005, 0.05534346867, 0.0580444660799995, 0.0610215936600004, 0.0643224600599991, 0.0679989058100032, 0.072122872869994, 0.0767768829000064, 0.0820773125999968, 0.0881712197999992)

msmc_model_expSFS_noDip=extrap_sso_model_trim_27_function(params_noDip,ns,pts_l)


ll_msmc_model_noDip = dadi.Inference.ll_multinom(msmc_model_expSFS_noDip, fs)
ll_msmc_model_noDip
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
def sso_model_trim_27_plusContraction((nu26, nu25, nu24, nu23, nu22, nu21, nu20, nu19, nu18, nu17, nu16, nu15, nu14, nu13, nu12, nu11, nu10, nu9, nu8, nu7, nu6, nu5, nu4, nu3, nu2, nu1, nu0, T26, T25, T24, T23, T22, T21, T20, T19, T18, T17, T16, T15, T14, T13, T12, T11, T10, T9, T8, T7, T6, T5, T4, T3, T2, T1, T0,contractionSize_RescByNanc,contractionGen_RescBy2Nanc),ns,pts):
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
 	### add contraction:
	phi = Integration.one_pop(phi, xx,contractionGen_RescBy2Nanc,  contractionSize_RescByNanc)
	# get expected SFS:
	fs = Spectrum.from_phi(phi,ns,(xx,))
	return fs


# wrap your function in Numerics.make_extrap_log_func to make extrapolation version of your function
extrap_sso_model_trim_27_plusContraction_function = Numerics.make_extrap_log_func(sso_model_trim_27_plusContraction)

params_plusContraction=(0.36675344664418, 0.141422067280315, 0.386870175412936, 0.609784494435127, 0.774150824662306, 0.91011303614748, 1.01928801965413, 1.09213414017289, 1.12514450221579, 1.12274330685933, 1.07023813979649, 1.07023813979649, 0.968579772770968, 0.968579772770968, 0.869604885707006, 0.869604885707006, 0.795226999398677, 0.795226999398677, 0.750983120146509, 0.750983120146509, 0.738213027247672, 0.738213027247672, 0.761765214135537, 0.761765214135537, 0.836551538729166, 0.836551538729166, 1, 0.030119771118, 0.030902139933, 0.0317262985630002, 0.0325958441059994, 0.0335133156900004, 0.0344855901200002, 0.0355149949300004, 0.0366078779399999, 0.0377705869699996, 0.03900946984, 0.0403308743700001, 0.0417474961999999, 0.0432646251800001, 0.0448991888299998, 0.04665859294, 0.0485650548800005, 0.0506312702899994, 0.0528826304500005, 0.05534346867, 0.0580444660799995, 0.0610215936600004, 0.0643224600599991, 0.0679989058100032, 0.072122872869994, 0.0767768829000064, 0.0820773125999968,0.0881712197999992,contractionSize_RescByNanc,contractionGen_RescBy2Nanc)

msmc_model_expSFS_wDip_wContraction=extrap_sso_model_trim_27_plusContraction_function(params_plusContraction,ns,pts_l)


ll_msmc_model_wDip_wContraction = dadi.Inference.ll_multinom(msmc_model_expSFS_wDip_wContraction, fs)

ll_msmc_model_wDip_wContraction

outputFigure=str(str(outdir)+"/"+str(pop)+"."+str(modelName)+".wDIP.AddRecentContraction.figure.png")
dadi.Plotting.plot_1d_comp_multinom(msmc_model_expSFS_wDip_wContraction, fs)
pyplot.title(str(modelName)+".wRecentContraction LL = "+str(ll_msmc_model_wDip_wContraction))
plt.savefig(outputFigure)

########################## try with No dip and dadi contraction ######
# converted that nu31 interval that was 0.0339730033926217 to  0.0929355795355233 (same as previous interval)
# converted that nu31 interval that was 0.141422067280315 to  0.386870175412936 (same as previous interval)
params_plusContraction_noDIP = (0.36675344664418, 0.386870175412936, 0.386870175412936, 0.609784494435127, 0.774150824662306, 0.91011303614748, 1.01928801965413, 1.09213414017289, 1.12514450221579, 1.12274330685933, 1.07023813979649, 1.07023813979649, 0.968579772770968, 0.968579772770968, 0.869604885707006, 0.869604885707006, 0.795226999398677, 0.795226999398677, 0.750983120146509, 0.750983120146509, 0.738213027247672, 0.738213027247672, 0.761765214135537, 0.761765214135537, 0.836551538729166, 0.836551538729166, 1, 0.030119771118, 0.030902139933, 0.0317262985630002, 0.0325958441059994, 0.0335133156900004, 0.0344855901200002, 0.0355149949300004, 0.0366078779399999, 0.0377705869699996, 0.03900946984, 0.0403308743700001, 0.0417474961999999, 0.0432646251800001, 0.0448991888299998, 0.04665859294, 0.0485650548800005, 0.0506312702899994, 0.0528826304500005, 0.05534346867, 0.0580444660799995, 0.0610215936600004, 0.0643224600599991, 0.0679989058100032, 0.072122872869994, 0.0767768829000064, 0.0820773125999968,0.0881712197999992,contractionSize_RescByNanc,contractionGen_RescBy2Nanc)

msmc_model_expSFS_wContraction_noDIP=extrap_sso_model_trim_27_plusContraction_function(params_plusContraction_noDIP,ns,pts_l)


ll_msmc_model_wContraction_noDIP = dadi.Inference.ll_multinom(msmc_model_expSFS_wContraction_noDIP, fs)

ll_msmc_model_wContraction_noDIP

outputFigure=str(str(outdir)+"/"+str(pop)+"."+str(modelName)+".noDIP.AddRecentContraction.figure.png")
dadi.Plotting.plot_1d_comp_multinom(msmc_model_expSFS_wContraction_noDIP, fs)
pyplot.title(str(modelName)+".noDIP.wRecentContraction LL = "+str(ll_msmc_model_wContraction_noDIP))
plt.savefig(outputFigure)

#################### Try to do dadi inference for that last epoch #######
# want to perturb initial parameters and do optimization

# Version of model with nu and T values plugged in plus one more epoch for dadi to optimize 
def sso_model_trim_27_plusContraction_forOptimization(params,ns,pts):
	nu,T=params
	xx = Numerics.default_grid(pts)
	# intialize phi with ancestral pop: (nu0 = 1)
	phi = PhiManip.phi_1D(xx,nu=1)
	# stays at nu0=1 for T0 duration of time:
	# followed by a number of time steps, with associated pop changes:
	phi = Integration.one_pop(phi, xx, 0.0881712197999992, 1)
	phi = Integration.one_pop(phi, xx, 0.0820773125999968, 0.836551538729166)
	phi = Integration.one_pop(phi, xx, 0.0767768829000064, 0.836551538729166)
	phi = Integration.one_pop(phi, xx, 0.072122872869994, 0.761765214135537)
	phi = Integration.one_pop(phi, xx, 0.0679989058100032, 0.761765214135537)
	phi = Integration.one_pop(phi, xx, 0.0643224600599991, 0.738213027247672)
	phi = Integration.one_pop(phi, xx, 0.0610215936600004, 0.738213027247672)
	phi = Integration.one_pop(phi, xx, 0.0580444660799995, 0.750983120146509)
	phi = Integration.one_pop(phi, xx, 0.05534346867, 0.750983120146509)
	phi = Integration.one_pop(phi, xx, 0.0528826304500005, 0.795226999398677)
	phi = Integration.one_pop(phi, xx, 0.0506312702899994, 0.795226999398677)
	phi = Integration.one_pop(phi, xx, 0.0485650548800005, 0.869604885707006)
	phi = Integration.one_pop(phi, xx, 0.04665859294, 0.869604885707006)
	phi = Integration.one_pop(phi, xx, 0.0448991888299998, 0.968579772770968)
	phi = Integration.one_pop(phi, xx, 0.0432646251800001, 0.968579772770968)
	phi = Integration.one_pop(phi, xx, 0.0417474961999999, 1.07023813979649)
	phi = Integration.one_pop(phi, xx, 0.0403308743700001, 1.07023813979649)
	phi = Integration.one_pop(phi, xx, 0.03900946984, 1.12274330685933)
	phi = Integration.one_pop(phi, xx, 0.0377705869699996, 1.12514450221579)
	phi = Integration.one_pop(phi, xx, 0.0366078779399999, 1.09213414017289)
	phi = Integration.one_pop(phi, xx, 0.0355149949300004, 1.01928801965413)
	phi = Integration.one_pop(phi, xx, 0.0344855901200002, 0.91011303614748)
	phi = Integration.one_pop(phi, xx, 0.0335133156900004, 0.774150824662306)
	phi = Integration.one_pop(phi, xx, 0.0325958441059994, 0.609784494435127)
	phi = Integration.one_pop(phi, xx, 0.0317262985630002, 0.386870175412936)
	phi = Integration.one_pop(phi, xx, 0.030902139933, 0.141422067280315)
	phi = Integration.one_pop(phi, xx, 0.030119771118, 0.36675344664418)
	### add contraction:
	phi = Integration.one_pop(phi, xx,T,  nu)
	# get expected SFS:
	fs = Spectrum.from_phi(phi,ns,(xx,))
	return fs


# wrap your function in Numerics.make_extrap_log_func to make extrapolation version of your function
extrap_sso_model_trim_27_plusContraction_forOptimization_function = Numerics.make_extrap_log_func(sso_model_trim_27_plusContraction_forOptimization)

upper_bound = [10, 2]
lower_bound = [1e-4, 1e-4]
p0 = [0.01,0.1] # initial parameters

# perturb parameters
### do 50 runs: 
p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound,                                  lower_bound=lower_bound)
# optimize:
print('Beginning optimization ************************************************')
popt = dadi.Inference.optimize_log(p0, fs, extrap_sso_model_trim_27_plusContraction_forOptimization_function, pts_l,
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=len(p0), maxiter=100)
print('Finshed optimization **************************************************')

# Calculate the best-fit model AFS.
model = extrap_sso_model_trim_27_plusContraction_forOptimization_function(popt, ns, pts_l)

# Likelihood of the data given the model AFS.
ll_model = dadi.Inference.ll_multinom(model, fs)
# calculate best fit theta
# calculate best fit theta
theta = dadi.Inference.optimal_sfs_scaling(model, fs)

###### model specific scaling of parameters (will depend on mu and L that you supply) #######
# param_names=("nuB","nuF","TF")
Nanc=theta / (4*mu*L)
nuB_scaled_dip=popt[0]*Nanc
nuF_scaled_dip=popt[1]*Nanc
#TF_scaled_gen=popt[2]*2*Nanc
scaled_param_names=("Nanc_FromTheta_scaled_dip","nuB_scaled_dip","nuF_scaled_dip")
scaled_popt=(Nanc,nuB_scaled_dip,nuF_scaled_dip)


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

outputSFS=str(outdir)+"/"+str(pop)+".dadi.inference."+str(modelName)+".runNum."+str(runNum)+".expSFS"

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


################### version without 

