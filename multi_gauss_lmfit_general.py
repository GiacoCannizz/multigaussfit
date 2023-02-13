'''
7 03 2019 Last version
Giacomo Cannizzaro
Multi-gaussian fit with lmfit package
https://lmfit.github.io/lmfit-py/index.html

To use: >python gaussfit.py file.dat
ascii file with 3 columns: WL, flux, flux_error

Always update list of prefixes, the printing of the results won't work otherwise  
'''

#*************************************#
#  REMEMBER TO UPDATE PREFIX LIST!!   #
#*************************************#

import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from lmfit import Model, Parameters
from scipy.constants import c

# Open file and normalise the flux
# file is opened from the command line

infile = open(sys.argv[1],"r")
xspec,flux,eflux=np.loadtxt(infile,unpack=True)

# flux normalisation - check your units and the plot label
# select WL interval
flux*=1.e17
eflux*=1.e17
xmin=4200.
xmax=5150
x=xspec[(xspec>xmin)&(xspec<xmax)]
y=flux[(xspec>xmin)&(xspec<xmax)]
ey=eflux[(xspec>xmin)&(xspec<xmax)]



# General definition of functions
# wl = x0 (centroid)
# x = x coordinate (wavelength for the spectrum)
def gauss(x,a,wl,sigma):
    return a*np.exp(-(x-wl)**2/(2*sigma**2))

def poly(x,slope,interc):
	return x*slope+interc


# Create a model with multiple gaussians and a 1st order polynomial
# to add a gaussian, you have to add Model(gauss, prefix='string')
# every gaussian will have a prefix (pref_) --> parameters are pref_a, pref_wl, pref_sigma
# calling 'make_params' will create parameters with the right prefix
#
# remember to update the list of prefixes, this is needed after for printing/plotting purposes
# you can add whatever component/model you want, there are pre-defined ones on the website
# or you can define your own and then use Model(f_name), using the prefix if you have more than one


mod = Model(poly) + \
	  Model(gauss, prefix='line1_') + \
	  Model(gauss, prefix='line2_') + \
	  Model(gauss, prefix='line3_')

params = mod.make_params()

#**********************************************************************#
# this list is important to update to get the right prints of results! #
#**********************************************************************#

prefixes = ['line1_','line2_','line3_']

# assigning initial values
# check https://lmfit.github.io/lmfit-py/constraints.html for example of possible constraints
# you can put min, max and set the parameter to vary or not, and tie them with mathematical expressions
# for every component you add to the model, you have to specify the components

params['slope'].set(value=-0.005)
params['interc'].set(value=4.)

# line1
params['line1_a'].set(value=1.)
params['line1_wl'].set(value=4334.,min=4300,max=4400, vary=True)
params['line1_sigma'].set(value=10.)
# line2
params['line2_a'].set(value=1.)
params['line2_wl'].set(value=4363.2)
params['line2_sigma'].set(value=5.)
# line3
params['line3_a'].set(value=1.)
params['line3_wl'].set(value=4958.9) 
params['line3_sigma'].set(value=5.)

# Examples of expressions to fix values
#  1) fixing wav separation between two lines (1 and 2)
#     params['line2_wl'].set(value=4363.2,expr='line1_wl-643.7')
#  2) Having the FWHM of two lines (1 and 2) to be the same
#     params['line2_sigma'].set(value=5., expr='line1_sigma/line1_wl*line2_wl')

# mod.fit does the fit
# prints a report of the fit, chisq ecc
# more options are available for the fit output, check the webpage 

result = mod.fit(y, params, x=x, weights=1/ey)
print result.fit_report(show_correl=False)



#*********************#
#  PRINTING RESULTS   #
#*********************#

#c=299792.458 #c in km/s - should use astropy constant
print '#line\t WV \t \t dWV \t FWHM(km/s) \t dFWHM(km/s) \t EW \t dEW \t flux \t \t dflux'

# values and errors for every parameter
# plus calculation of EW and propagation
for pref in prefixes:
	slop=result.params['slope'].value
	dslop=result.params['slope'].stderr
	dinter=result.params['interc'].stderr
	WL=result.params[pref+'wl'].value
	dWL=result.params[pref+'wl'].stderr
	sig=abs(result.params[pref+'sigma'].value)
	dsig=result.params[pref+'sigma'].stderr
	cont=result.params['slope'].value*WL+result.params['interc'].value
	FWHM=2.3548*sig
	dFWHM=2.3548*dsig
	FWkms=FWHM*c/WL
	dFWkms=c*2.3548*np.sqrt( (dsig/WL)**2 + (sig*dWL/WL**2)**2 ) 
	amp=result.params[pref+'a'].value
	damp=result.params[pref+'a'].stderr
	flux=2.50663*amp*sig
	dflux=flux*np.sqrt( (damp/amp)**2 + (dsig/sig)**2 )
	EW=flux/cont
	dEW=EW*np.sqrt( (damp/amp)**2 + (dsig/sig)**2 + (dinter/cont)**2 + (dWL*slop/cont)**2 + (WL*dslop/cont)**2 )
	print pref,'\t{0} \t {4} \t {1} \t {5} \t {2} \t {6} \t {3} \t {7} \\\\'\
		.format('%.2f'%WL, '%.2f'%FWkms, '%.2f'%EW, '%.3f'%flux, '%.2f'%dWL, '%.2f'%dFWkms, '%.2f'%dEW, '%.2f'%dflux)	
	

#*************#
#  PLOTTING   #
#*************#

# to have the single components for plot
comps = result.eval_components()

# 'auto scale' the plot depending on the flux min and max 
plmin=min(y)
plmax=max(y)

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=17)


fig=plt.figure(figsize=(8, 6))
gs1=gridspec.GridSpec(20,1)
ax1 = fig.add_subplot(gs1[:17,:])
ax2 = fig.add_subplot(gs1[17:,:])
#ax1.set_ylim(4.5,15.)
ax1.set_ylim(plmin-0.6,plmax+2)
ax2.set_ylim(-2.5,2.5)
ax1.tick_params(labelbottom='off')
ax1.set_ylabel('Flux [1E-17 erg/$\mathrm{cm^2}$/s/\AA]')
ax2.set_xlabel('Restframe Wavelength [\AA]')
ax1.plot(x,y,'k',lw=1,zorder=100)
ax1.plot(x,result.best_fit,'m',lw=2)

# single components
for pref in prefixes:
	ax1.plot(x, comps[pref]+comps['poly'],linestyle='--')
	
#best fit	
ax2.plot(x,result.best_fit-y,'k', linestyle=':')

# this if you want to how the 3-sigma 'shade' around the best model fit
# dely = result.eval_uncertainty(sigma=3)
# ax1.fill_between(x, result.best_fit-dely, result.best_fit+dely, color="#ABABAB")

title = infile.name
plt.savefig("fit_"+title+".eps", dpi=200)
plt.show()

