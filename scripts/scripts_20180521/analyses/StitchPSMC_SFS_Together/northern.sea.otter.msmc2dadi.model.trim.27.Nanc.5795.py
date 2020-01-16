# dadi function for nso msmc model; trim point = 27; Nanc=5795
# note: order of variables goes from oldest --> youngest; so nu0 is Nanc, T0 is time for Nanc. nu32 would be most recent if there were 33 time indices.
# This is reverse order for time index numbering in MSMC so don't get confused. They are listed in descending order because that's how they come from MSMC.

from dadi import Numerics, PhiManip, Integration, Spectrum, Plotting, Inference
import numpy as np
import matplotlib.pyplot as pyplot


Nanc = 5795
# params in dadi units (nu = Nx/Nanc; T = gen/(2*Nanc)) in same order that they should go into function (all nus from nu32 --> nu0 followed by Ts from T32 --> T0) -- note that order is from recent --> old because that is how it comes out of msmc. But then dadi function starts with oldest (nu0) first.
params = (0.66278005763842, 0.486537327256331, 0.418811341101633, 0.454087486635882, 0.57893158551086, 0.751518688070356, 0.923258438398994, 1.05264451217583, 1.12120902352939, 1.13323637139827, 1.07810101738927, 1.07810101738927, 0.96238333076211, 0.96238333076211, 0.853951454898572, 0.853951454898572, 0.776471091610656, 0.776471091610656, 0.73277072167976, 0.73277072167976, 0.722873076226689, 0.722873076226689, 0.751275450679311, 0.751275450679311, 0.831977694051747, 0.831977694051747, 1, 0.0325326988008, 0.0333777959977, 0.0342675085774996, 0.0352072263180003, 0.0361993446760004, 0.037247357026, 0.0383602463309993, 0.0395410069120003, 0.0407966255179997, 0.0421340888980007, 0.0435623800149997, 0.0450914799390003, 0.04673136974, 0.0484950248090002, 0.0503974167509998, 0.052455513385, 0.0546882787439995, 0.0571186692890003, 0.0597776263369998, 0.0626930968840031, 0.0659109938519982, 0.0694782282700026, 0.0734407130599937, 0.0779022513500069, 0.0829327106299926, 0.0886518637400001, 0.0952293888700031)



def nso_model_trim_27((nu26, nu25, nu24, nu23, nu22, nu21, nu20, nu19, nu18, nu17, nu16, nu15, nu14, nu13, nu12, nu11, nu10, nu9, nu8, nu7, nu6, nu5, nu4, nu3, nu2, nu1, nu0, T26, T25, T24, T23, T22, T21, T20, T19, T18, T17, T16, T15, T14, T13, T12, T11, T10, T9, T8, T7, T6, T5, T4, T3, T2, T1, T0),ns,pts):
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
extrap_nso_model_trim_27_function = Numerics.make_extrap_log_func(nso_model_trim_27)


# Version of model with nu and T values plugged in plus one more epoch for dadi to optimize 
def nso_model_trim_27_plusContraction_forOptimization((contractionGen_RescBy2Nanc,contractionSize_RescByNanc),ns,pts):
	xx = Numerics.default_grid(pts)
	# intialize phi with ancestral pop: (nu0 = 1)
	phi = PhiManip.phi_1D(xx,nu=1)
	# stays at nu0=1 for T0 duration of time:
	# followed by a number of time steps, with associated pop changes:
	phi = Integration.one_pop(phi, xx, 0.0952293888700031, 1)
	phi = Integration.one_pop(phi, xx, 0.0886518637400001, 0.831977694051747)
	phi = Integration.one_pop(phi, xx, 0.0829327106299926, 0.831977694051747)
	phi = Integration.one_pop(phi, xx, 0.0779022513500069, 0.751275450679311)
	phi = Integration.one_pop(phi, xx, 0.0734407130599937, 0.751275450679311)
	phi = Integration.one_pop(phi, xx, 0.0694782282700026, 0.722873076226689)
	phi = Integration.one_pop(phi, xx, 0.0659109938519982, 0.722873076226689)
	phi = Integration.one_pop(phi, xx, 0.0626930968840031, 0.73277072167976)
	phi = Integration.one_pop(phi, xx, 0.0597776263369998, 0.73277072167976)
	phi = Integration.one_pop(phi, xx, 0.0571186692890003, 0.776471091610656)
	phi = Integration.one_pop(phi, xx, 0.0546882787439995, 0.776471091610656)
	phi = Integration.one_pop(phi, xx, 0.052455513385, 0.853951454898572)
	phi = Integration.one_pop(phi, xx, 0.0503974167509998, 0.853951454898572)
	phi = Integration.one_pop(phi, xx, 0.0484950248090002, 0.96238333076211)
	phi = Integration.one_pop(phi, xx, 0.04673136974, 0.96238333076211)
	phi = Integration.one_pop(phi, xx, 0.0450914799390003, 1.07810101738927)
	phi = Integration.one_pop(phi, xx, 0.0435623800149997, 1.07810101738927)
	phi = Integration.one_pop(phi, xx, 0.0421340888980007, 1.13323637139827)
	phi = Integration.one_pop(phi, xx, 0.0407966255179997, 1.12120902352939)
	phi = Integration.one_pop(phi, xx, 0.0395410069120003, 1.05264451217583)
	phi = Integration.one_pop(phi, xx, 0.0383602463309993, 0.923258438398994)
	phi = Integration.one_pop(phi, xx, 0.037247357026, 0.751518688070356)
	phi = Integration.one_pop(phi, xx, 0.0361993446760004, 0.57893158551086)
	phi = Integration.one_pop(phi, xx, 0.0352072263180003, 0.454087486635882)
	phi = Integration.one_pop(phi, xx, 0.0342675085774996, 0.418811341101633)
	phi = Integration.one_pop(phi, xx, 0.0333777959977, 0.486537327256331)
	phi = Integration.one_pop(phi, xx, 0.0325326988008, 0.66278005763842)
	### add contraction:
	phi = Integration.one_pop(phi, xx,contractionGen_RescBy2Nanc,  contractionSize_RescByNanc)
	# get expected SFS:
	fs = Spectrum.from_phi(phi,ns,(xx,))
	return fs


# wrap your function in Numerics.make_extrap_log_func to make extrapolation version of your function
extrap_nso_model_trim_27_plusContraction_forOptimization_function = Numerics.make_extrap_log_func(nso_model_trim_27_plusContraction_forOptimization)
