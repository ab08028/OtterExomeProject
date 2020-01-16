# dadi function for nso msmc model; trim point = 31; Nanc=12244
# note: order of variables goes from oldest --> youngest; so nu0 is Nanc, T0 is time for Nanc. nu32 would be most recent if there were 33 time indices.
# This is reverse order for time index numbering in MSMC so don't get confused. They are listed in descending order because that's how they come from MSMC.

from dadi import Numerics, PhiManip, Integration, Spectrum, Plotting, Inference
import numpy as np
import matplotlib.pyplot as pyplot


Nanc = 12244
# params in dadi units (nu = Nx/Nanc; T = gen/(2*Nanc)) in same order that they should go into function (all nus from nu32 --> nu0 followed by Ts from T32 --> T0) -- note that order is from recent --> old because that is how it comes out of msmc. But then dadi function starts with oldest (nu0) first.
params = (0.313698420919823, 0.230281508201515, 0.198226326897981, 0.214922772457405, 0.274012354630085, 0.355699033219889, 0.436984654092703, 0.498223984644426, 0.53067604576918, 0.536368673347274, 0.51027272599619, 0.51027272599619, 0.455502738352359, 0.455502738352359, 0.404181175725737, 0.404181175725737, 0.367509179736121, 0.367509179736121, 0.346825490052124, 0.346825490052124, 0.342140865471663, 0.342140865471663, 0.35558390726732, 0.35558390726732, 0.393780841557747, 0.393780841557747, 0.473306970094387, 0.473306970094387, 0.636426341832739, 0.636426341832739, 1, 0.0153979530984, 0.0157979434920999, 0.0162190506574998, 0.0166638256140001, 0.0171334021480001, 0.0176294336979999, 0.0181561719629996, 0.0187150341760001, 0.0193093272139998, 0.0199423579540003, 0.0206183780949998, 0.021342111747, 0.0221182830199999, 0.022953033257, 0.0238534486229998, 0.0248275601049999, 0.0258843435119997, 0.027034664297, 0.0282931672009998, 0.0296730797320013, 0.031196132795999, 0.0328845297100011, 0.0347600013799969, 0.0368716785500031, 0.0392526299899963, 0.0419595450199999, 0.0450727335100013, 0.0486819535500004, 0.0529194802199982, 0.0579648297000023, 0.0640778280399997)



def nso_model_trim_31((nu30, nu29, nu28, nu27, nu26, nu25, nu24, nu23, nu22, nu21, nu20, nu19, nu18, nu17, nu16, nu15, nu14, nu13, nu12, nu11, nu10, nu9, nu8, nu7, nu6, nu5, nu4, nu3, nu2, nu1, nu0, T30, T29, T28, T27, T26, T25, T24, T23, T22, T21, T20, T19, T18, T17, T16, T15, T14, T13, T12, T11, T10, T9, T8, T7, T6, T5, T4, T3, T2, T1, T0),ns,pts):
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
	# get expected SFS:
	fs = Spectrum.from_phi(phi,ns,(xx,))
	return fs


# wrap your function in Numerics.make_extrap_log_func to make extrapolation version of your function
extrap_nso_model_trim_31_function = Numerics.make_extrap_log_func(nso_model_trim_31)


# Version of model with nu and T values plugged in plus one more epoch for dadi to optimize 
def nso_model_trim_31_plusContraction_forOptimization((contractionGen_RescBy2Nanc,contractionSize_RescByNanc),ns,pts):
	xx = Numerics.default_grid(pts)
	# intialize phi with ancestral pop: (nu0 = 1)
	phi = PhiManip.phi_1D(xx,nu=1)
	# stays at nu0=1 for T0 duration of time:
	# followed by a number of time steps, with associated pop changes:
	phi = Integration.one_pop(phi, xx, 0.0640778280399997, 1)
	phi = Integration.one_pop(phi, xx, 0.0579648297000023, 0.636426341832739)
	phi = Integration.one_pop(phi, xx, 0.0529194802199982, 0.636426341832739)
	phi = Integration.one_pop(phi, xx, 0.0486819535500004, 0.473306970094387)
	phi = Integration.one_pop(phi, xx, 0.0450727335100013, 0.473306970094387)
	phi = Integration.one_pop(phi, xx, 0.0419595450199999, 0.393780841557747)
	phi = Integration.one_pop(phi, xx, 0.0392526299899963, 0.393780841557747)
	phi = Integration.one_pop(phi, xx, 0.0368716785500031, 0.35558390726732)
	phi = Integration.one_pop(phi, xx, 0.0347600013799969, 0.35558390726732)
	phi = Integration.one_pop(phi, xx, 0.0328845297100011, 0.342140865471663)
	phi = Integration.one_pop(phi, xx, 0.031196132795999, 0.342140865471663)
	phi = Integration.one_pop(phi, xx, 0.0296730797320013, 0.346825490052124)
	phi = Integration.one_pop(phi, xx, 0.0282931672009998, 0.346825490052124)
	phi = Integration.one_pop(phi, xx, 0.027034664297, 0.367509179736121)
	phi = Integration.one_pop(phi, xx, 0.0258843435119997, 0.367509179736121)
	phi = Integration.one_pop(phi, xx, 0.0248275601049999, 0.404181175725737)
	phi = Integration.one_pop(phi, xx, 0.0238534486229998, 0.404181175725737)
	phi = Integration.one_pop(phi, xx, 0.022953033257, 0.455502738352359)
	phi = Integration.one_pop(phi, xx, 0.0221182830199999, 0.455502738352359)
	phi = Integration.one_pop(phi, xx, 0.021342111747, 0.51027272599619)
	phi = Integration.one_pop(phi, xx, 0.0206183780949998, 0.51027272599619)
	phi = Integration.one_pop(phi, xx, 0.0199423579540003, 0.536368673347274)
	phi = Integration.one_pop(phi, xx, 0.0193093272139998, 0.53067604576918)
	phi = Integration.one_pop(phi, xx, 0.0187150341760001, 0.498223984644426)
	phi = Integration.one_pop(phi, xx, 0.0181561719629996, 0.436984654092703)
	phi = Integration.one_pop(phi, xx, 0.0176294336979999, 0.355699033219889)
	phi = Integration.one_pop(phi, xx, 0.0171334021480001, 0.274012354630085)
	phi = Integration.one_pop(phi, xx, 0.0166638256140001, 0.214922772457405)
	phi = Integration.one_pop(phi, xx, 0.0162190506574998, 0.198226326897981)
	phi = Integration.one_pop(phi, xx, 0.0157979434920999, 0.230281508201515)
	phi = Integration.one_pop(phi, xx, 0.0153979530984, 0.313698420919823)
	### add contraction:
	phi = Integration.one_pop(phi, xx,contractionGen_RescBy2Nanc,  contractionSize_RescByNanc)
	# get expected SFS:
	fs = Spectrum.from_phi(phi,ns,(xx,))
	return fs


# wrap your function in Numerics.make_extrap_log_func to make extrapolation version of your function
extrap_nso_model_trim_31_plusContraction_forOptimization_function = Numerics.make_extrap_log_func(nso_model_trim_31_plusContraction_forOptimization)
