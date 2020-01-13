# dadi function for sso msmc model; trim point = 27; Nanc=5467
# note: order of variables goes from oldest --> youngest; so nu0 is Nanc, T0 is time for Nanc. nu32 would be most recent if there were 33 time indices.
# This is reverse order for time index numbering in MSMC so don't get confused. They are listed in descending order because that's how they come from MSMC.

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
	phi = Integration.one_pop(phi,xx,T0,nu0)
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


# Version of model with nu and T values plugged in plus one more epoch for dadi to optimize 
def sso_model_trim_27_plusContraction_forOptimization((contractionGen_RescBy2Nanc,contractionSize_RescByNanc),ns,pts):
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
	phi = Integration.one_pop(phi, xx,contractionGen_RescBy2Nanc,  contractionSize_RescByNanc)
	# get expected SFS:
	fs = Spectrum.from_phi(phi,ns,(xx,))
	return fs


# wrap your function in Numerics.make_extrap_log_func to make extrapolation version of your function
extrap_sso_model_trim_27_plusContraction_forOptimization_function = Numerics.make_extrap_log_func(sso_model_trim_27_plusContraction_forOptimization)
