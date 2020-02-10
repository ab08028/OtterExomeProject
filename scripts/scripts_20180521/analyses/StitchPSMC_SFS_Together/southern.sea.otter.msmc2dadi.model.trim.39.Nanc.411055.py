# dadi function for sso msmc model; trim point = 39; Nanc=411055
# note: order of variables goes from oldest --> youngest; so nu0 is Nanc, T0 is time for Nanc. nu32 would be most recent if there were 33 time indices.
# This is reverse order for time index numbering in MSMC so don't get confused. They are listed in descending order because that's how they come from MSMC.

from dadi import Numerics, PhiManip, Integration, Spectrum, Plotting, Inference
import numpy as np
import matplotlib.pyplot as pyplot


Nanc = 411055
# params in dadi units (nu = Nx/Nanc; T = gen/(2*Nanc)) in same order that they should go into function (all nus from nu32 --> nu0 followed by Ts from T32 --> T0) -- note that order is from recent --> old because that is how it comes out of msmc. But then dadi function starts with oldest (nu0) first.
params = (0.00487809781986972, 0.00188102029958802, 0.00514566550504812, 0.00811059429737348, 0.0102967906221188, 0.0121051907162397, 0.0135573004479985, 0.014526208865738, 0.014965271610991, 0.0149333338993196, 0.0142349755244366, 0.0142349755244366, 0.0128828424685752, 0.0128828424685752, 0.01156640172282, 0.01156640172282, 0.0105771196632592, 0.0105771196632592, 0.00998864265534717, 0.00998864265534717, 0.00981879077556433, 0.00981879077556433, 0.0101320526482532, 0.0101320526482532, 0.0111267672454692, 0.0111267672454692, 0.0133007552199023, 0.0133007552199023, 0.0179070175840453, 0.0179070175840453, 0.0283778878420236, 0.0283778878420236, 0.0553680897108008, 0.0553680897108008, 0.137643054169846, 0.137643054169846, 0.447459632029814, 0.447459632029814, 1, 0.000400615702920001, 0.00041102179902, 0.000421983731220002, 0.000433549343639992, 0.000445752408600005, 0.000458684392800002, 0.000472376254200005, 0.000486912423599999, 0.000502377331799994, 0.0005188554096, 0.000536431087800002, 0.000555273227999998, 0.000575452189200002, 0.000597193120199998, 0.0006205945236, 0.000645951907200007, 0.000673434132599993, 0.000703378923000007, 0.0007361099298, 0.000772035235199994, 0.000811633280400006, 0.000855537296399988, 0.000904436801400043, 0.00095928867779992, 0.00102119052600009, 0.00109169024399996, 0.00117274381199999, 0.00126646200000007, 0.00137692562999991, 0.00150807480600003, 0.00166722686399998, 0.00186366919200004, 0.00211288077000004, 0.00243920581199994, 0.00288500043600006, 0.00353089605600002, 0.00455208658199998, 0.00641589649199996, 0.0109679830740001)



def sso_model_trim_39((nu38, nu37, nu36, nu35, nu34, nu33, nu32, nu31, nu30, nu29, nu28, nu27, nu26, nu25, nu24, nu23, nu22, nu21, nu20, nu19, nu18, nu17, nu16, nu15, nu14, nu13, nu12, nu11, nu10, nu9, nu8, nu7, nu6, nu5, nu4, nu3, nu2, nu1, nu0, T38, T37, T36, T35, T34, T33, T32, T31, T30, T29, T28, T27, T26, T25, T24, T23, T22, T21, T20, T19, T18, T17, T16, T15, T14, T13, T12, T11, T10, T9, T8, T7, T6, T5, T4, T3, T2, T1, T0),ns,pts):
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
extrap_sso_model_trim_39_function = Numerics.make_extrap_log_func(sso_model_trim_39)


# Version of model with nu and T values plugged in plus one more epoch for dadi to optimize 
def sso_model_trim_39_plusContraction_forOptimization((contractionGen_RescBy2Nanc,contractionSize_RescByNanc),ns,pts):
	xx = Numerics.default_grid(pts)
	# intialize phi with ancestral pop: (nu0 = 1)
	phi = PhiManip.phi_1D(xx,nu=1)
	# stays at nu0=1 for T0 duration of time:
	# followed by a number of time steps, with associated pop changes:
	phi = Integration.one_pop(phi, xx, 0.0109679830740001, 1)
	phi = Integration.one_pop(phi, xx, 0.00641589649199996, 0.447459632029814)
	phi = Integration.one_pop(phi, xx, 0.00455208658199998, 0.447459632029814)
	phi = Integration.one_pop(phi, xx, 0.00353089605600002, 0.137643054169846)
	phi = Integration.one_pop(phi, xx, 0.00288500043600006, 0.137643054169846)
	phi = Integration.one_pop(phi, xx, 0.00243920581199994, 0.0553680897108008)
	phi = Integration.one_pop(phi, xx, 0.00211288077000004, 0.0553680897108008)
	phi = Integration.one_pop(phi, xx, 0.00186366919200004, 0.0283778878420236)
	phi = Integration.one_pop(phi, xx, 0.00166722686399998, 0.0283778878420236)
	phi = Integration.one_pop(phi, xx, 0.00150807480600003, 0.0179070175840453)
	phi = Integration.one_pop(phi, xx, 0.00137692562999991, 0.0179070175840453)
	phi = Integration.one_pop(phi, xx, 0.00126646200000007, 0.0133007552199023)
	phi = Integration.one_pop(phi, xx, 0.00117274381199999, 0.0133007552199023)
	phi = Integration.one_pop(phi, xx, 0.00109169024399996, 0.0111267672454692)
	phi = Integration.one_pop(phi, xx, 0.00102119052600009, 0.0111267672454692)
	phi = Integration.one_pop(phi, xx, 0.00095928867779992, 0.0101320526482532)
	phi = Integration.one_pop(phi, xx, 0.000904436801400043, 0.0101320526482532)
	phi = Integration.one_pop(phi, xx, 0.000855537296399988, 0.00981879077556433)
	phi = Integration.one_pop(phi, xx, 0.000811633280400006, 0.00981879077556433)
	phi = Integration.one_pop(phi, xx, 0.000772035235199994, 0.00998864265534717)
	phi = Integration.one_pop(phi, xx, 0.0007361099298, 0.00998864265534717)
	phi = Integration.one_pop(phi, xx, 0.000703378923000007, 0.0105771196632592)
	phi = Integration.one_pop(phi, xx, 0.000673434132599993, 0.0105771196632592)
	phi = Integration.one_pop(phi, xx, 0.000645951907200007, 0.01156640172282)
	phi = Integration.one_pop(phi, xx, 0.0006205945236, 0.01156640172282)
	phi = Integration.one_pop(phi, xx, 0.000597193120199998, 0.0128828424685752)
	phi = Integration.one_pop(phi, xx, 0.000575452189200002, 0.0128828424685752)
	phi = Integration.one_pop(phi, xx, 0.000555273227999998, 0.0142349755244366)
	phi = Integration.one_pop(phi, xx, 0.000536431087800002, 0.0142349755244366)
	phi = Integration.one_pop(phi, xx, 0.0005188554096, 0.0149333338993196)
	phi = Integration.one_pop(phi, xx, 0.000502377331799994, 0.014965271610991)
	phi = Integration.one_pop(phi, xx, 0.000486912423599999, 0.014526208865738)
	phi = Integration.one_pop(phi, xx, 0.000472376254200005, 0.0135573004479985)
	phi = Integration.one_pop(phi, xx, 0.000458684392800002, 0.0121051907162397)
	phi = Integration.one_pop(phi, xx, 0.000445752408600005, 0.0102967906221188)
	phi = Integration.one_pop(phi, xx, 0.000433549343639992, 0.00811059429737348)
	phi = Integration.one_pop(phi, xx, 0.000421983731220002, 0.00514566550504812)
	phi = Integration.one_pop(phi, xx, 0.00041102179902, 0.00188102029958802)
	phi = Integration.one_pop(phi, xx, 0.000400615702920001, 0.00487809781986972)
	### add contraction:
	phi = Integration.one_pop(phi, xx,contractionGen_RescBy2Nanc,  contractionSize_RescByNanc)
	# get expected SFS:
	fs = Spectrum.from_phi(phi,ns,(xx,))
	return fs


# wrap your function in Numerics.make_extrap_log_func to make extrapolation version of your function
extrap_sso_model_trim_39_plusContraction_forOptimization_function = Numerics.make_extrap_log_func(sso_model_trim_39_plusContraction_forOptimization)
