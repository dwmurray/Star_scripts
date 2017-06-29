import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
from matplotlib import rcParams
import scipy.integrate as si

def shu_derivs( x, y) :
	v = y[0]
	alpha = y[1] 
	dvdx = (alpha*(x-v) - 2e0/x)*(x-v)/((x-v)**2. - 1.)
	dalphadx = alpha*(alpha-2./x*(x-v))*(x-v)/((x-v)**2. - 1.)
	return [dvdx, dalphadx]

def shu_solution( x0, x1, v0, alpha0 ) : 
	integrator = si.ode( shu_derivs).set_integrator( "dopri5")
	integrator.set_initial_value( [v0,alpha0], x0)
	dx = 0.001*(x1-x0)
	vs = []
	alphas = [] 
	xs = []
	while( integrator.successful() and integrator.t > x1) : 
		result = integrator.integrate(integrator.t+dx)
		xs.append(integrator.t)
		vs.append( result[0])
		alphas.append( result[1])

	xs = np.array(xs)
	vs = np.array(vs)
	alphas = np.array(alphas)
	return xs, vs, alphas

import argparse
parser = argparse.ArgumentParser(description = "start number to end number, stride length, reduction method")
parser.add_argument('--fit', action='store_true')
parser.add_argument('--mdot', action='store_true')
parser.add_argument('--density', action='store_true')
args = parser.parse_args()

mp = 1.6e-24
Msun = 2e33
# rhobin is in number density
if args.density:
    value = 'density'
if args.mdot:
    value = 'mdot'
files = ['avg_{0}_1.0Msun_shellsphere_36_sinks.txt'.format(value), 'avg_{0}_4.0Msun_shellsphere_30_sinks.txt'.format(value)]
labels = ["$1.0 M_\odot$", "$4.0 M_\odot$"]
ltypes = ["solid", "solid"]
lweights = [2, 2]

plt.rc("text", usetex=True)
plt.rc('font', family='serif')
rcParams.update({'figure.autolayout': True})
rcParams['axes.linewidth']    = 2
rcParams['lines.linewidth']   = 2
rcParams['xtick.major.width'] = 3    # major tick width in points
rcParams['xtick.major.size']  = 14    # major tick width in points
rcParams['xtick.minor.size']  = 8    # major tick width in points
rcParams['xtick.minor.width'] = 2    # major tick width in points
rcParams['ytick.major.width'] = 3    # major tick width in points
rcParams['ytick.minor.width'] = 2    # major tick width in points
rcParams['ytick.major.size']  = 14    # major tick width in points
rcParams['ytick.minor.size']  = 8    # major tick width in points
rcParams['xtick.labelsize']   = 25


for file, label, ltype, lw in zip(files,labels,ltypes,lweights): 
    rbin, avgbin = np.loadtxt(file, usecols=[0,1], unpack=True, comments="=")

    lr = np.log10(rbin)
    lavg = np.log10(avgbin)

    lr_min = np.log10(5.0e-3)
    if args.fit:
        if True:
            lr_max = np.log10(4.0e-1)
            lr_mask = lr[lr>lr_min]
            lavg_mask = lavg[lr>lr_min]
            p = np.polyfit( lr_mask[lr_mask<lr_max], lavg_mask[lr_mask<lr_max], 1)
        #	k = np.polyfit( lt[lt>=lt_kink], lm[lt>=lt_kink], 1)
            plt.loglog(rbin, avgbin, linewidth=lw,label="{0} $\\alpha_1={1:2.3f}$".format(label, -1.0*p[0]),ls=ltype)
        #	pl.loglog( tarray, mtot, linewidth=lw,label="{0} $\\alpha_1={1:2.2f}$ $\\alpha_2={2:2.2f}$".format(label, -1.0*p[0], -1.0*k[0]),ls=ltype)
            lt =  np.arange(lr_min, lr_max, 0.1 )
            plt.loglog( 1e1**lt, 1e1**(p[0]*lt + p[1]), ls="dotted", linewidth=1,color="black")
        elif False:
            p = np.polyfit( lr[lr>lr_min], lavg[lr>lr_min], 1)
            plt.loglog(rbin, avgbin, linewidth=lw,label="{0} $\\alpha_1={1:2.3f}$".format(label, -1.0*p[0]),ls=ltype)
    #    plt.loglog(rbin, avgbin, linewidth=lw,label="{0} $\\alpha_1={1}$".format(label, '7'),ls=ltype)
    else:
        plt.loglog(rbin, avgbin, linewidth=lw,label="{0}".format(label),ls=ltype)

if args.density:
    plt.ylabel('$N$ $({\\rm cm^{-3}})$', fontsize=25)
if args.mdot:
    plt.ylabel('$\dot{M}$ $({\\rm M_\odot \, yr^{-1}})$', fontsize=25)
    for A in np.arange( 2.001, 20., 1.) :
	x0 = 100.
	alpha0 = A/x0**2 - A*(A-2.)/(2.*x0**4) 
	v0 = -(A-2.)/x0 - ( 1. - A/6.)*(A-2.)/x0**3
	cs = 2.65e4
	x, v, alpha = shu_solution( x0, 0.01, v0, alpha0) 
	mdot = -cs**3./6.67e-8*x*x* alpha* v/2e33*3.15e7
	if A == 4.001:
            plt.loglog( x, mdot, label="A={0}".format(A), linestyle='-' )


plt.xlabel('$r$ $({\\rm pc})$', fontsize=25)
plt.xlim(3e-3, 4e0)	
plt.legend(loc='best')
plt.rc('text', usetex=False)
plt.savefig('avg_overplot.pdf')
