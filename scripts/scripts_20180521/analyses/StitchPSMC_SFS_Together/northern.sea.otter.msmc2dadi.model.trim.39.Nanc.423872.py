# dadi function for nso msmc model; trim point = 39; Nanc=423872
# note: order of variables goes from oldest --> youngest; so nu0 is Nanc, T0 is time for Nanc. nu32 would be most recent if there were 33 time indices.
# This is reverse order for time index numbering in MSMC so don't get confused. They are listed in descending order because that's how they come from MSMC.

from dadi import Numerics, PhiManip, Integration, Spectrum, Plotting, Inference
import numpy as np
import matplotlib.pyplot as pyplot


Nanc = 423872
# params in dadi units (nu = Nx/Nanc; T = gen/(2*Nanc)) in same order that they should go into function (all nus from nu32 --> nu0 followed by Ts from T32 --> T0) -- note that order is from recent --> old because that is how it comes out of msmc. But then dadi function starts with oldest (nu0) first.
params = (0.00906164920249148, 0.00665202661532087, 0.00572606464444715, 0.00620836650667638, 0.00791525767814158, 0.0102748998584465, 0.0126229568853081, 0.014391946761725, 0.015329373201259, 0.015493813283558, 0.0147399927198178, 0.0147399927198178, 0.0131578795124961, 0.0131578795124961, 0.0116753792318683, 0.0116753792318683, 0.0106160536470003, 0.0106160536470003, 0.0100185742603333, 0.0100185742603333, 0.00988325185587543, 0.00988325185587543, 0.0102715742727033, 0.0102715742727033, 0.0113749499866631, 0.0113749499866631, 0.0136721814394649, 0.0136721814394649, 0.018384129049815, 0.018384129049815, 0.0288864992559446, 0.0288864992559446, 0.0557012298410961, 0.0557012298410961, 0.137194521186969, 0.137194521186969, 0.444441187582237, 0.444441187582237, 1, 0.00044479296072, 0.00045634728293, 0.000468511594749995, 0.000481359586200005, 0.000494924008400006, 0.0005092526234, 0.00052446824789999, 0.000540611820800003, 0.000557778866199996, 0.00057606490820001, 0.000595592763499996, 0.000616498895100004, 0.000638919765999999, 0.000663032778100003, 0.000689042625899997, 0.000717181296499999, 0.000747708069599993, 0.000780936810100004, 0.000817290553299997, 0.000857151395600042, 0.000901147066799975, 0.000949918943000035, 0.00100409475399991, 0.00106509371500009, 0.0011338710669999, 0.001212064366, 0.00130199348300004, 0.00140625121500002, 0.00152865852599995, 0.00167440101000007, 0.001850984132, 0.00206932493199995, 0.00234593543300002, 0.00270810823500001, 0.00320305953600002, 0.00392030906399995, 0.00505418013100008, 0.00712323213699992, 0.012177548731)



def nso_model_trim_39((nu38, nu37, nu36, nu35, nu34, nu33, nu32, nu31, nu30, nu29, nu28, nu27, nu26, nu25, nu24, nu23, nu22, nu21, nu20, nu19, nu18, nu17, nu16, nu15, nu14, nu13, nu12, nu11, nu10, nu9, nu8, nu7, nu6, nu5, nu4, nu3, nu2, nu1, nu0, T38, T37, T36, T35, T34, T33, T32, T31, T30, T29, T28, T27, T26, T25, T24, T23, T22, T21, T20, T19, T18, T17, T16, T15, T14, T13, T12, T11, T10, T9, T8, T7, T6, T5, T4, T3, T2, T1, T0),ns,pts):
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
	phi = Integration.one_pop(phi, xx, T33, nu33)
	phi = Integration.one_pop(phi, xx, T34, nu34)
	phi = Integration.one_pop(phi, xx, T35, nu35)
	phi = Integration.one_pop(phi, xx, T36, nu36)
	phi = Integration.one_pop(phi, xx, T37, nu37)
	phi = Integration.one_pop(phi, xx, T38, nu38)
	# get expected SFS:
	fs = Spectrum.from_phi(phi,ns,(xx,))
	return fs


# wrap your function in Numerics.make_extrap_log_func to make extrapolation version of your function
extrap_nso_model_trim_39_function = Numerics.make_extrap_log_func(nso_model_trim_39)


# Version of model with nu and T values plugged in plus one more epoch for dadi to optimize 
def nso_model_trim_39_plusContraction_forOptimization((contractionGen_RescBy2Nanc,contractionSize_RescByNanc),ns,pts):
	xx = Numerics.default_grid(pts)
	# intialize phi with ancestral pop: (nu0 = 1)
	phi = PhiManip.phi_1D(xx,nu=1)
	# stays at nu0=1 for T0 duration of time:
	# followed by a number of time steps, with associated pop changes:
	phi = Integration.one_pop(phi, xx, 0.012177548731, 1)
	phi = Integration.one_pop(phi, xx, 0.00712323213699992, 0.444441187582237)
	phi = Integration.one_pop(phi, xx, 0.00505418013100008, 0.444441187582237)
	phi = Integration.one_pop(phi, xx, 0.00392030906399995, 0.137194521186969)
	phi = Integration.one_pop(phi, xx, 0.00320305953600002, 0.137194521186969)
	phi = Integration.one_pop(phi, xx, 0.00270810823500001, 0.0557012298410961)
	phi = Integration.one_pop(phi, xx, 0.00234593543300002, 0.0557012298410961)
	phi = Integration.one_pop(phi, xx, 0.00206932493199995, 0.0288864992559446)
	phi = Integration.one_pop(phi, xx, 0.001850984132, 0.0288864992559446)
	phi = Integration.one_pop(phi, xx, 0.00167440101000007, 0.018384129049815)
	phi = Integration.one_pop(phi, xx, 0.00152865852599995, 0.018384129049815)
	phi = Integration.one_pop(phi, xx, 0.00140625121500002, 0.0136721814394649)
	phi = Integration.one_pop(phi, xx, 0.00130199348300004, 0.0136721814394649)
	phi = Integration.one_pop(phi, xx, 0.001212064366, 0.0113749499866631)
	phi = Integration.one_pop(phi, xx, 0.0011338710669999, 0.0113749499866631)
	phi = Integration.one_pop(phi, xx, 0.00106509371500009, 0.0102715742727033)
	phi = Integration.one_pop(phi, xx, 0.00100409475399991, 0.0102715742727033)
	phi = Integration.one_pop(phi, xx, 0.000949918943000035, 0.00988325185587543)
	phi = Integration.one_pop(phi, xx, 0.000901147066799975, 0.00988325185587543)
	phi = Integration.one_pop(phi, xx, 0.000857151395600042, 0.0100185742603333)
	phi = Integration.one_pop(phi, xx, 0.000817290553299997, 0.0100185742603333)
	phi = Integration.one_pop(phi, xx, 0.000780936810100004, 0.0106160536470003)
	phi = Integration.one_pop(phi, xx, 0.000747708069599993, 0.0106160536470003)
	phi = Integration.one_pop(phi, xx, 0.000717181296499999, 0.0116753792318683)
	phi = Integration.one_pop(phi, xx, 0.000689042625899997, 0.0116753792318683)
	phi = Integration.one_pop(phi, xx, 0.000663032778100003, 0.0131578795124961)
	phi = Integration.one_pop(phi, xx, 0.000638919765999999, 0.0131578795124961)
	phi = Integration.one_pop(phi, xx, 0.000616498895100004, 0.0147399927198178)
	phi = Integration.one_pop(phi, xx, 0.000595592763499996, 0.0147399927198178)
	phi = Integration.one_pop(phi, xx, 0.00057606490820001, 0.015493813283558)
	phi = Integration.one_pop(phi, xx, 0.000557778866199996, 0.015329373201259)
	phi = Integration.one_pop(phi, xx, 0.000540611820800003, 0.014391946761725)
	phi = Integration.one_pop(phi, xx, 0.00052446824789999, 0.0126229568853081)
	phi = Integration.one_pop(phi, xx, 0.0005092526234, 0.0102748998584465)
	phi = Integration.one_pop(phi, xx, 0.000494924008400006, 0.00791525767814158)
	phi = Integration.one_pop(phi, xx, 0.000481359586200005, 0.00620836650667638)
	phi = Integration.one_pop(phi, xx, 0.000468511594749995, 0.00572606464444715)
	phi = Integration.one_pop(phi, xx, 0.00045634728293, 0.00665202661532087)
	phi = Integration.one_pop(phi, xx, 0.00044479296072, 0.00906164920249148)
	### add contraction:
	phi = Integration.one_pop(phi, xx,contractionGen_RescBy2Nanc,  contractionSize_RescByNanc)
	# get expected SFS:
	fs = Spectrum.from_phi(phi,ns,(xx,))
	return fs


# wrap your function in Numerics.make_extrap_log_func to make extrapolation version of your function
extrap_nso_model_trim_39_plusContraction_forOptimization_function = Numerics.make_extrap_log_func(nso_model_trim_39_plusContraction_forOptimization)
