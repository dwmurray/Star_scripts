import numpy as np
import matplotlib.pyplot as plt
import matplotlib
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
parser = argparse.ArgumentParser(description = "What do you want to plot?")
parser.add_argument('--fit', action='store_true')
parser.add_argument('--mdot', action='store_true')
parser.add_argument('--density', action='store_true')
parser.add_argument('--shu', action='store_true')
parser.add_argument('--kink', action='store_true')
parser.add_argument('--nocone', action='store_true')
args = parser.parse_args()

mp = 1.6e-24
Msun = 2e33

if args.kink :
	args.fit = True

# rhobin is in number density
if args.density:
    value = 'density'
    files = ['jet/avg_{0}_1.0Msun_shellsphere_36_sinks.txt'.format(value), 'jet/avg_{0}_4.0Msun_shellsphere_30_sinks.txt'.format(value), \
		     'jet/avg_density_by_hand_shellsphere_80kyr.txt', \
		     'nojet/avg_{0}_1.0Msun_shellsphere_9_sinks.txt'.format(value), 'nojet/avg_{0}_4.0Msun_shellsphere_23_sinks.txt'.format(value)]
    labels = ["$1.0 M_\odot$ jet", "$4.0 M_\odot$ jet","-80 kyr jet", "$1.0 M_\odot$ no jet", "$4.0 M_\odot$ no jet"]
    ltypes = ["solid", "solid", ":", "--", "--"]
    lweights = [2, 4, 3, 2, 4]
    lcolors = ['b', 'b', 'b', 'r', 'r']
    if args.nocone:
	    files = ['jet/avg_{0}_1.0Msun_shellsphere_36_sinks.txt'.format(value), 'jet/avg_{0}_4.0Msun_shellsphere_30_sinks.txt'.format(value), \
			     'jet/nojetremoval/avg_{0}_1.0Msun_shellsphere_28_sinks.txt'.format(value), 'jet/nojetremoval/avg_{0}_4.0Msun_shellsphere_17_sinks.txt'.format(value)]
	    labels = ["$1.0 M_\odot$ jet", "$4.0 M_\odot$ jet", "$1.0 M_\odot$ no jet sub", "$4.0 M_\odot$ no jet sub"]
	    ltypes = ["solid", "solid", "--", "--"]
	    lweights = [2, 4, 2, 4]
	    lcolors = ['b', 'b', 'r', 'r']

if args.mdot:
    value = 'mdot'
    files = ['jet/avg_{0}_1.0Msun_shellsphere_36_sinks.txt'.format(value), 'jet/avg_{0}_4.0Msun_shellsphere_30_sinks.txt'.format(value), \
		     'nojet/avg_{0}_1.0Msun_shellsphere_9_sinks.txt'.format(value), 'nojet/avg_{0}_4.0Msun_shellsphere_23_sinks.txt'.format(value), \
		     'jet/avg_mdot_by_hand_shellsphere_80kyr.txt']
    labels = ["$Jet \, 1.0 M_\odot$", "$Jet \, 4.0 M_\odot$", \
		      "$No \, Jet \, 1.0 M_\odot$", "$No \, Jet \, 4.0 M_\odot$", \
		      "$Jet \, -80 kyr$"]
    ltypes = ["solid", "solid", "dashed", "dashed", ":"]
    lweights = [2, 4, 2, 4, 2]
    lcolors = ['b', 'b', 'r', 'r', 'b']


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
rcParams['ytick.labelsize']   = 25


if args.density:
	rcParams['axes.linewidth']    = 1
	rcParams['lines.linewidth']   = 1
	rcParams['xtick.minor.width'] = 1    # major tick width in points
	rcParams['ytick.minor.width'] = 1    # major tick width in points

	rcParams['xtick.major.width'] = 2    # major tick width in points
	rcParams['ytick.major.width'] = 2    # major tick width in points

	rcParams['xtick.major.size']  = 10    # major tick width in points
	rcParams['ytick.major.size']  = 10    # major tick width in points

	rcParams['xtick.minor.size']  = 4   # major tick width in points
	rcParams['ytick.minor.size']  = 4    # major tick width in points
	rcParams['xtick.labelsize']   = 20
	rcParams['ytick.labelsize']   = 20

	#plt.ylabel('$N$ $({\\rm cm^{-3}})$', fontsize=25)
	host = host_subplot(111, axes_class=AA.Axes)
	par1 = host.twinx()
#	host.set_xlim(3e-3, 3e0)
	Ndensity_Y_min = 1e1
	Ndensity_Y_max = 1e8
	#host.set_ylim(Ndensity_Y_min, Ndensity_Y_max)

	Mdensity_Y_min = Ndensity_Y_min * 2. * mp
	Mdensity_Y_max = Ndensity_Y_max * 2. * mp
	#par1.set_ylim(Mdensity_Y_min, Mdensity_Y_max)
	par1.set_yscale('log')
	par1.set_ylabel('$n$ $({\\rm cm^{-3}})$', fontsize=20)
	host.set_xlabel('$r$ $({\\rm pc})$', fontsize=20)
	host.set_ylabel('$\\rho$ $({\\rm g \\, cm}^{-3})$', fontsize=20)

	host.axis["left"].label.set_fontsize(20)
	host.axis["bottom"].label.set_fontsize(20)
	par1.axis["right"].label.set_fontsize(20)

for file, label, ltype, lw, lc in zip(files,labels,ltypes,lweights, lcolors): 
    rbin, avgbin = np.loadtxt(file, usecols=[0,1], unpack=True, comments="=")

    lr = np.log10(rbin)
    lavg = np.log10(avgbin)

    lr_min = np.log10(2.0e-2)
    if args.fit:
	    if label != "-80 kyr jet":
		    lr_max = np.log10(1.0e0)
		    lr_mask = lr[lr>lr_min]
		    lavg_mask = lavg[lr>lr_min]
		    p = np.polyfit( lr_mask[lr_mask<lr_max], lavg_mask[lr_mask<lr_max], 1)
		    plt.loglog(rbin, avgbin, linewidth=lw,label="{0} $\\alpha_1={1:2.3f}$".format(label, -1.0*p[0]),ls=ltype, color=lc)
#	    if args.kink:
		    #k = np.polyfit( lr[lr>=lr_ma], lm[lt>=lt_kink], 1)
		    #pl.loglog( tarray, mtot, linewidth=lw,label="{0} $\\alpha_1={1:2.2f}$ $\\alpha_2={2:2.2f}$".format(label, -1.0*p[0], -1.0*k[0]),ls=ltype)

		    #If you want to plot the fits.
		    #lt =  np.arange(lr_min, lr_max, 0.1 )
		    #plt.loglog( 1e1**lt, 1e1**(p[0]*lt + p[1]), ls="dotted", linewidth=1,color="black")
	    else:
		    plt.loglog(rbin, avgbin, linewidth=lw,label="{0}".format(label),ls=ltype, color=lc)
    elif args.mdot:
	    plt.loglog(rbin, avgbin, color=lc, linewidth=lw,label="{0}".format(label),ls=ltype)

if args.mdot:
    plt.ylabel('$\dot{M}$ $({\\rm M_\odot \, yr^{-1}})$', fontsize=25)
    plt.ylim(1e-6, 1e-3)
    if args.shu:
	    for A in np.arange( 2.001, 7., 0.5) :
		    x0 = 10.
		    alpha0 = A/x0**2 - A*(A-2.)/(2.*x0**4) 
		    v0 = -(A-2.)/x0 - ( 1. - A/6.)*(A-2.)/x0**3
		    cs = 2.65e4
		    x, v, alpha = shu_solution( x0, 0.001, v0, alpha0) 
		    mdot = -cs**3./6.67e-8*x*x* alpha* v/2e33*3.15e7
		    if A == 3.501:
			    plt.loglog( x, mdot, color='g', ls='dashdot', label="Shu".format(A))

plt.xlabel('$r$ $({\\rm pc})$', fontsize=25)
plt.xlim(3e-3, 4e0)	
if args.density:
	plt.legend(loc='best')
plt.rc('text', usetex=False)
plt.savefig('avg_{0}.pdf'.format(value))
