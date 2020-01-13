# dadi function for sso msmc model; trim point = 33; Nanc=22759
# note: order of variables goes from oldest --> youngest; so nu0 is Nanc, T0 is time for Nanc. nu32 would be most recent if there were 33 time indices.
# This is reverse order for time index numbering in MSMC so don't get confused. They are listed in descending order because that's how they come from MSMC.

from dadi import Numerics, PhiManip, Integration, Spectrum, Plotting, Inference
import numpy as np
import matplotlib.pyplot as pyplot


Nanc = 22759
# params in dadi units (nu = Nx/Nanc; T = gen/(2*Nanc)) in same order that they should go into function (all nus from nu32 --> nu0 followed by Ts from T32 --> T0) -- note that order is from recent --> old because that is how it comes out of msmc. But then dadi function starts with oldest (nu0) first.
params = (0.0881030544009925, 0.0339730033926217, 0.0929355795355233, 0.146484994149822, 0.185969764821238, 0.218631178707224, 0.244857652102702, 0.262357053342665, 0.270286941253667, 0.269710116013023, 0.257097104104347, 0.257097104104347, 0.232676303911964, 0.232676303911964, 0.208900140554491, 0.208900140554491, 0.191032772098617, 0.191032772098617, 0.18040432146964, 0.18040432146964, 0.177336636081359, 0.177336636081359, 0.182994441404337, 0.182994441404337, 0.200959926621754, 0.200959926621754, 0.240224202954715, 0.240224202954715, 0.323417652253806, 0.323417652253806, 0.512531459731541, 0.512531459731541, 1, 0.00723549801, 0.00742344193499999, 0.00762142478500003, 0.00783031066999984, 0.00805070955000007, 0.00828427340000002, 0.00853156135000007, 0.00879409829999996, 0.00907340914999987, 0.00937101879999997, 0.00968845215000001, 0.0100287589999999, 0.0103932101, 0.0107858718499999, 0.0112085233, 0.0116665016000001, 0.0121628565499998, 0.0127036877500001, 0.01329484065, 0.0139436855999999, 0.0146588637000001, 0.0154518116999998, 0.0163349829500007, 0.0173256596499985, 0.0184436655000015, 0.0197169569999992, 0.0211808609999998, 0.0228735000000012, 0.0248685774999983, 0.0272372555000004, 0.0301116919999995, 0.0336596260000006, 0.0381606225000007)



def sso_model_trim_33((nu32, nu31, nu30, nu29, nu28, nu27, nu26, nu25, nu24, nu23, nu22, nu21, nu20, nu19, nu18, nu17, nu16, nu15, nu14, nu13, nu12, nu11, nu10, nu9, nu8, nu7, nu6, nu5, nu4, nu3, nu2, nu1, nu0, T32, T31, T30, T29, T28, T27, T26, T25, T24, T23, T22, T21, T20, T19, T18, T17, T16, T15, T14, T13, T12, T11, T10, T9, T8, T7, T6, T5, T4, T3, T2, T1, T0),ns,pts):
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


# Version of model with nu and T values plugged in plus one more epoch for dadi to optimize 
def sso_model_trim_33_plusContraction_forOptimization((contractionGen_RescBy2Nanc,contractionSize_RescByNanc),ns,pts):
	xx = Numerics.default_grid(pts)
	# intialize phi with ancestral pop: (nu0 = 1)
	phi = PhiManip.phi_1D(xx,nu=1)
	# stays at nu0=1 for T0 duration of time:
	# followed by a number of time steps, with associated pop changes:
	phi = Integration.one_pop(phi, xx, 0.0381606225000007, 1)
	phi = Integration.one_pop(phi, xx, 0.0336596260000006, 0.512531459731541)
	phi = Integration.one_pop(phi, xx, 0.0301116919999995, 0.512531459731541)
	phi = Integration.one_pop(phi, xx, 0.0272372555000004, 0.323417652253806)
	phi = Integration.one_pop(phi, xx, 0.0248685774999983, 0.323417652253806)
	phi = Integration.one_pop(phi, xx, 0.0228735000000012, 0.240224202954715)
	phi = Integration.one_pop(phi, xx, 0.0211808609999998, 0.240224202954715)
	phi = Integration.one_pop(phi, xx, 0.0197169569999992, 0.200959926621754)
	phi = Integration.one_pop(phi, xx, 0.0184436655000015, 0.200959926621754)
	phi = Integration.one_pop(phi, xx, 0.0173256596499985, 0.182994441404337)
	phi = Integration.one_pop(phi, xx, 0.0163349829500007, 0.182994441404337)
	phi = Integration.one_pop(phi, xx, 0.0154518116999998, 0.177336636081359)
	phi = Integration.one_pop(phi, xx, 0.0146588637000001, 0.177336636081359)
	phi = Integration.one_pop(phi, xx, 0.0139436855999999, 0.18040432146964)
	phi = Integration.one_pop(phi, xx, 0.01329484065, 0.18040432146964)
	phi = Integration.one_pop(phi, xx, 0.0127036877500001, 0.191032772098617)
	phi = Integration.one_pop(phi, xx, 0.0121628565499998, 0.191032772098617)
	phi = Integration.one_pop(phi, xx, 0.0116665016000001, 0.208900140554491)
	phi = Integration.one_pop(phi, xx, 0.0112085233, 0.208900140554491)
	phi = Integration.one_pop(phi, xx, 0.0107858718499999, 0.232676303911964)
	phi = Integration.one_pop(phi, xx, 0.0103932101, 0.232676303911964)
	phi = Integration.one_pop(phi, xx, 0.0100287589999999, 0.257097104104347)
	phi = Integration.one_pop(phi, xx, 0.00968845215000001, 0.257097104104347)
	phi = Integration.one_pop(phi, xx, 0.00937101879999997, 0.269710116013023)
	phi = Integration.one_pop(phi, xx, 0.00907340914999987, 0.270286941253667)
	phi = Integration.one_pop(phi, xx, 0.00879409829999996, 0.262357053342665)
	phi = Integration.one_pop(phi, xx, 0.00853156135000007, 0.244857652102702)
	phi = Integration.one_pop(phi, xx, 0.00828427340000002, 0.218631178707224)
	phi = Integration.one_pop(phi, xx, 0.00805070955000007, 0.185969764821238)
	phi = Integration.one_pop(phi, xx, 0.00783031066999984, 0.146484994149822)
	phi = Integration.one_pop(phi, xx, 0.00762142478500003, 0.0929355795355233)
	phi = Integration.one_pop(phi, xx, 0.00742344193499999, 0.0339730033926217)
	phi = Integration.one_pop(phi, xx, 0.00723549801, 0.0881030544009925)
	### add contraction:
	phi = Integration.one_pop(phi, xx,contractionGen_RescBy2Nanc,  contractionSize_RescByNanc)
	# get expected SFS:
	fs = Spectrum.from_phi(phi,ns,(xx,))
	return fs


# wrap your function in Numerics.make_extrap_log_func to make extrapolation version of your function
extrap_sso_model_trim_33_plusContraction_forOptimization_function = Numerics.make_extrap_log_func(sso_model_trim_33_plusContraction_forOptimization)
