import matplotlib
matplotlib.use("Agg")
from matplotlib import rc
rc("text", usetex=True)
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
rcParams['axes.linewidth'] = 2
rcParams['lines.linewidth'] = 2
rcParams['xtick.major.width'] = 3    # major tick width in points
rcParams['xtick.major.size'] = 14    # major tick width in points
rcParams['xtick.minor.size'] = 8    # major tick width in points
rcParams['xtick.minor.width'] = 2    # major tick width in points
rcParams['ytick.major.width'] = 3    # major tick width in points
rcParams['ytick.minor.width'] = 2    # major tick width in points
rcParams['ytick.major.size'] = 14    # major tick width in points
rcParams['ytick.minor.size'] = 8    # major tick width in points
rcParams['xtick.labelsize']=25
import matplotlib.backends.backend_pdf as mat_pdf
import argparse
import glob
import sys
import os
import math
import h5py
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA



import numpy as np
import scipy.integrate as si
import matplotlib.pyplot as pl 


mp = 1.6e-24
pi = np.pi
parsec = 3e18
sec_in_yr = np.pi* 1e7
Msun = 2e33
G = 6.67e-8
c_s = 2.65e4
cs = c_s
Text_64_filename = 'ramses_jeans_32/rad_profile_0171_shellsphere_1.0.out'
Text_32_filename = 'ramses_jeans_32c/rad_profile_0150_shellsphere_1.0.out'
Text_16_filename = 'ramses_jeans_32c/rad_profile_0145_shellsphere_1.out'
Text_08_filename = 'ramses_jeans_32c/rad_profile_0131_shellsphere_1.out'
Text_04_filename = 'ramses_jeans_32c/rad_profile_0105_shellsphere_1.out'


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

pl.clf()
for A in np.arange( 2.001, 10.001, 1.) :
	x0 = 1000.
	alpha0 = A/x0**2 - A*(A-2.)/(2.*x0**4) 
	v0 = -(A-2.)/x0 - ( 1. - A/6.)*(A-2.)/x0**3

	x, v, alpha = shu_solution( x0, 0.01, v0, alpha0) 
	mdot = -cs**3./6.67e-8*x*x* alpha* v/2e33*3.15e7
#	print type(A)
	if A > 4.501 and A < 5.501:
#		print 'plotting'
		#pl.loglog( x/120., mdot, '--', label="A={0}".format(A))
		pl.loglog( x/120., mdot, '--', label="Shu")



rbin64, vrbin64, vrmsbin64, vrmsnbin64, vKbin64, vmagbin64, vmagnbin64, mTbin64, rhobin64, mdotbin64, norm64, angXbin64, angYbin64, angZbin64, vphi_magbin64, sum_of_velocities64, partMass64, part_creation_time64, current_time64, vrms_r_bin64, vrms_l_bin64, vrms_theta_bin64, vrms_phi_bin64 = np.loadtxt(Text_64_filename, unpack=True)
rbin4, vrbin4, vrmsbin4, vrmsnbin4, vKbin4, vmagbin4, vmagnbin4, mTbin4, rhobin4, mdotbin4, norm4, angXbin4, angYbin4, angZbin4, vphi_magbin4, sum_of_velocities4, partMass4, part_creation_time4, current_time4, vrms_r_bin4, vrms_l_bin4, vrms_theta_bin4, vrms_phi_bin4 = np.loadtxt(Text_04_filename, unpack=True)
rbin8, vrbin8, vrmsbin8, vrmsnbin8, vKbin8, vmagbin8, vmagnbin8, mTbin8, rhobin8, mdotbin8, norm8, angXbin8, angYbin8, angZbin8, vphi_magbin8, sum_of_velocities8, partMass8, part_creation_time8, current_time8, vrms_r_bin8, vrms_l_bin8, vrms_theta_bin8, vrms_phi_bin8 = np.loadtxt(Text_08_filename, unpack=True)
rbin16, vrbin16, vrmsbin16, vrmsnbin16, vKbin16, vmagbin16, vmagnbin16, mTbin16, rhobin16, mdotbin16, norm16, angXbin16, angYbin16, angZbin16, vphi_magbin16, sum_of_velocities16, partMass16, part_creation_time16, current_time16, vrms_r_bin16, vrms_l_bin16, vrms_theta_bin16, vrms_phi_bin16 = np.loadtxt(Text_16_filename, unpack=True)
rbin32, vrbin32, vrmsbin32, vrmsnbin32, vKbin32, vmagbin32, vmagnbin32, mTbin32, rhobin32, mdotbin32, norm32, angXbin32, angYbin32, angZbin32, vphi_magbin32, sum_of_velocities32, partMass32, part_creation_time32, current_time32, vrms_r_bin32, vrms_l_bin32, vrms_theta_bin32, vrms_phi_bin32 = np.loadtxt(Text_32_filename, unpack=True)


pl.rc('text', usetex=False)
pl.rc('font', family='serif')
filestar = 146
toutkyrs = 6.146
pl.loglog( rbin4, -4.0*pi*rhobin4*vrbin4*(parsec*rbin4)**2*sec_in_yr/Msun, label='$(t-t_*)={0:4.1f}$'.format((105-filestar)*toutkyrs))
pl.loglog( rbin8, -4.0*pi*rhobin8*vrbin8*(parsec*rbin8)**2*sec_in_yr/Msun, '.', label='$(t-t_*)={0:4.1f}$'.format((131-filestar)*toutkyrs))
pl.loglog( rbin16, -4.0*pi*rhobin16*vrbin16*(parsec*rbin16)**2*sec_in_yr/Msun, '-^', label='$(t-t_*)={0:4.1f}$'.format((145-filestar)*toutkyrs))
pl.loglog( rbin32, -4.0*pi*rhobin32*vrbin32*(parsec*rbin32)**2*sec_in_yr/Msun, '.-', label='$(t-t_*)={0:4.1f}$'.format((150-filestar)*toutkyrs))# label='$v_T150$')
pl.loglog( rbin64, -4.0*pi*rhobin64*vrbin64*(parsec*rbin64)**2*sec_in_yr/Msun, label='$(t-t_*)={0:4.1f}$'.format((171-filestar)*toutkyrs)) #label='$v_T171$')

pl.gcf().subplots_adjust(bottom=0.15)
pl.ylabel(r'$\dot{M}$ $({\rm M_{\odot}\, yr}^{-1})$', fontsize=25)
pl.xlabel(r'$r$ $({\rm pc})$', fontsize=28)
pl.xlim( 2e-3, 3.0) 
pl.ylim( 1.e-7, 2.e-3) 
pl.xticks(fontsize=24)
pl.yticks(fontsize=24)
pl.gcf().subplots_adjust(bottom=0.15)
#pl.rc('text', usetex=False)
#pl.legend(loc="best") 
pl.legend(loc="best", fontsize=18, frameon=False, ncol=2)
pl.savefig('Ramses_mdot_shu.pdf')
