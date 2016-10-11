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
#xtick.minor.width    : 0.5    # minor tick width in points

import matplotlib.pyplot as pl
import numpy as np
import matplotlib.backends.backend_pdf as mat_pdf
import argparse
import glob
import sys
import os
import math
import h5py

from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA


#def plotting_routine(rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, rhobin, mdotbin, norm, angXbin, angYbin, angZbin, vphi_magbin):
#	velocity_plot(rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, vphi_magbin)
#	pdf_timestep.savefig()
#	density_plot(rbin, mTbin, rhobin)
#	pdf_timestep.savefig()
#	ang_moment_plot(rbin, mTbin, vphi_magbin)
#	pdf_timestep.savefig()
#	total_mass_plot(rbin, mTbin)
#	pdf_timestep.savefig()
#	Mdot_plot(rbin, vrbin, rhobin, mdotbin)
#	pdf_timestep.savefig()
#	pdf_close(pdf_timestep)
#
def velocity_plot(rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, vphi_magbin):

	index_vk_min=vKbin.argmin()
	index_turb_min=vrmsbin[20:65].argmin()
	index_turb_min = index_turb_min + 20
	MC15_test_velocity = np.sqrt(vrmsbin*vrmsbin + vphi_magbin*vphi_magbin)

	# use this to figure out what % of vkep vr infall is.
	vr_vk_ratio = vrbin / vKbin
	vr_vk_ratio = abs(vr_vk_ratio)
	print vr_vk_ratio
	print 'Minimum % vr is of vk: ', vr_vk_ratio.min()

	# Now plot
	pl.clf()
	pl.rc('text', usetex=True)
	pl.rc('font', family='serif')
	pl.loglog( rbin, -vrbin/1e5, 'bv-', label='$u_r$')
	pl.loglog( rbin, vrmsbin/1e5, 'g.-', label='$v_T$')  
	#pl.loglog( rbin, vmagbin/1e5, 'm.', label='$v_{Total}$')  
	pl.loglog( rbin, vKbin/1e5, 'r--', label='$v_K$')  
	#pl.loglog( rbin, MC15_test_velocity/1e5, '.--', color='orange', label='v_t + v_phi')  
	pl.loglog( rbin, vphi_magbin/1e5, 'k+', label='$v_\phi$')
	#This is the sum
	#pl.loglog( rbin, np.sqrt(sum_of_velocities**2 + c_s**2)/1e5, 'c', label='$\sqrt{v_r^2 + v_T^2 + v_\phi^2}$')  
	# Need to mask this so that if not fully refined still got values.
#	power_law_fit = least_squares_fitting(log_x_sorted=rbin, log_y_sorted=vrmsbin)	
#	pl.loglog( rbin, power_law_fit, 'k.-', label='$powerlaw$')
#	print power_law_fit
#	sys.exit()
	pl.axhline( c_s/1e5, color = 'black',label='$c_s$')  
#	if not (withNoParticle) or (withdisk):
#		if disk_radius_index != 0 or particleMass != 0.0:
#			vy_second_point = vrmsbin[index_turb_min] / 1e5
#			vx_plot = [r_star, r_star]
#			vy_plot = [1e-2, vy_second_point]
#			pl.plot (vx_plot, vy_plot, '.--', color = 'g')
#			pl.plot (r_star, vy_second_point, color = 'g', marker = 'o')
#			#pl.plot (rbin[index_turb_min], vrmsbin[index_turb_min] , color = 'g', marker = 'o')
#			pl.annotate('$R_*$ = ' + "%.2f"%r_star + 'pc', xy=(r_star, 1e-2), xytext=(r_star, 1.25e-2), fontsize=20)
#	pl.legend(loc=3, fontsize=21, frameon=False, ncol=2)
	pl.legend(loc='best', fontsize=21, frameon=False, ncol=2)
	pl.ylim(1e-2, 3e0)
	pl.xlim(3e-3, 3e0)

	pl.xticks(fontsize=24)
	pl.yticks(fontsize=24)
	#pl.minorticks_on()
	#pl.tick_params('both',length=8, width=1, which='minor')
	#pl.tick_params('both',length=10, width=1.5, which='major')
	pl.gcf().subplots_adjust(bottom=0.15)
	pl.ylabel(r'$u_r, \, v_T, \, v_K, \, v_\phi, \, c_s$ $({\rm km\, s}^{-1})$', fontsize = 25)
	pl.xlabel(r'$r$ $({\rm pc})$', fontsize=25)
	#pl.title('{0} framestep {1:04d}'.format(vel_Plot_label, framestep))
	pl.rc('text', usetex=False)
	#pl.title('Particle {0} {1:3d} kyr after formation'.format(particle_name, int(partAge)))

def Turbulent_velocity_plot(rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, vphi_magbin):

	#vrms_subtract_bin = vrmsbin - vrms_r_bin
	#vrms_total = np.sqrt(vrms_r_bin**2 + vrms_l_bin**2)
	#vrms_total = vrms_r_bin + vrms_l_bin
	# Now plot
	pl.clf()
	pl.rc('text', usetex=True)
	pl.rc('font', family='serif')
	#pl.loglog( rbin, -vrbin/1e5, 'bv-', label='$v_r$')
	pl.loglog( rbin, vrmsbin/1e5, 'go-', label='$v_T$', lw = 1)  
	pl.loglog( rbin, vrms_r_bin/1e5, '-', color='blue', label='$v_{T,r}$')  
	pl.loglog( rbin, vrms_l_bin/1e5, 'r--', label='$v_{T,l}$', lw = 2)  
	#pl.loglog( rbin, np.sqrt(2.0*(vrms_l_bin/1e5)**2 + (vrms_r_bin/1e5)**2), '^-', label='$v_{all}$', lw = 2)  

	#pl.loglog( rbin, vrms_theta_bin/1e5, '-', color = 'brown', label=r'$v_{\theta}$', lw = 3)  
	#pl.loglog( rbin, vrms_phi_bin/1e5, 'bv-', label=r'$v_{\phi}$')  
	#pl.loglog( rbin, vrms_subtract_bin/1e5, 'o-', color = 'grey', lw = 0.5, label=r'$v_{Tl2}$')  
	#pl.loglog( rbin, vrms_total/1e5, 'o-', color = 'grey', lw = 0.5, label=r'$v_{Total}$')  
	#pl.loglog( rbin, vKbin/1e5, 'r--', label='$v_K$')  
	#pl.loglog( rbin, vphi_magbin/1e5, 'k+', label='$v_\phi$')
	pl.axhline( c_s/1e5, color = 'black',label='$c_s$')  
	pl.legend(loc=0, fontsize=22, frameon=False)
	pl.ylim(1e-1, 3e0)
	pl.xlim(3e-3, 3e0)

	pl.xticks(fontsize=24)
	pl.yticks(fontsize=24)
	#pl.minorticks_on()
	#pl.tick_params('both',length=8, width=1, which='minor')
	#pl.tick_params('both',length=10, width=1.5, which='major')
	pl.gcf().subplots_adjust(bottom=0.15)
	pl.ylabel(' $v_T, \\, v_{T,r}, \\, v_{T,l}$ $({\\rm km\\, s}^{-1})$', fontsize=25)
	pl.xlabel('$r$ $(\\rm pc)$', fontsize=25)
	#pl.title('{0} framestep {1:04d}'.format(vel_Plot_label, framestep))
	pl.rc('text', usetex=False)
	#pl.title('Particle {0} {1:3d} kyr after formation'.format(particle_name, int(partAge)))
	
def density_plot(rbin, mTbin, rhobin):
	pl.clf()
	pl.rc('text', usetex=True)
	pl.rc('font', family='serif')
	for i in range(len(rbin)):
		print rbin[i], rhobin[i]
	host = host_subplot(111, axes_class=AA.Axes)
	par1 = host.twinx()
	Ndensity_Y_min = 1e1
	Ndensity_Y_max = 1e8
	host.set_xlim(3e-3, 4e0)
	host.set_ylim(Ndensity_Y_min, Ndensity_Y_max)
	Mdensity_Y_min = Ndensity_Y_min * 2. * mp
	Mdensity_Y_max = Ndensity_Y_max * 2. * mp
	par1.set_ylim(Mdensity_Y_min, Mdensity_Y_max)
	par1.set_yscale('log')
	#pl.xticks(fontsize=24)
	pl.xticks(fontsize=30)
	#pl.yticks(fontsize=24)
	pl.yticks(fontsize=30)
	#par1.yticks(fontsize= 20)
	#pl.minorticks_on()
	#pl.tick_params('both',length=8, width=1, which='minor')
	#pl.tick_params('both',length=10, width=1.5, which='major')
	#par1.tick_params('both',length=10, width=1.5, which='major')
	pl.gcf().subplots_adjust(bottom=0.15)
	host.set_ylabel('$n$ $({\\rm cm^{-3}})$', fontsize=25)
	host.set_xlabel('$r$ $({\\rm pc})$', fontsize=25)
	par1.set_ylabel('$\\rho$ $({\\rm g \\, cm}^{-3})$', fontsize=25)
        host.axis["left"].label.set_fontsize(25)
        host.axis["bottom"].label.set_fontsize(25)
        par1.axis["right"].label.set_fontsize(25)
	# rhobin should be in g/cm**3
	host.loglog(rbin, rhobin/mp, label='Number density') 	
	#host.loglog(rbin, dlnrho_dlnr, label='derivative')
	#pl.legend()
	# If plotting multiple things, this may be useful
	#host.axis["left"].label.set_color(p1.get_color())
	#par1.axis["right"].label.set_color(p2.get_color())
#	if not (withNoParticle) or (withdisk):
#		if disk_radius_index != 0 or particleMass != 0.0:
#			x_plot = [r_star, r_star]
#			y_plot = [1e0, 1e4]
#			host.plot (x_plot, y_plot, '.--', color = 'g')
#			host.plot (r_star, chosen_ratio_number * disk_particle_mass, color = 'g', marker = 'o')
#			host.annotate('$R_*$ = ' + "%.2f"%r_star + 'pc', xy=(r_star, 1e32), xytext=(r_star, 1.25e32))
#	pl.rc('text', usetex=True)
	pl.rc('text', usetex=False)
	#pl.title('Particle {0} {1:03d} kyr after formation'.format(particle_name, int(partAge)))
	#pl.title('Run of Density {0:03d} kyr after formation'.format(int(partAge/1e3)))
	
def ang_moment_plot(rbin, mTbin, vphi_magbin):
	pl.clf()
	pl.rc('text', usetex=True)
	pl.rc('font', family='serif')
	#pl.loglog( rbin, (rbin*parsec)*vphi_magbin/np.sqrt(G*mTbin*(parsec*rbin)))
	pl.loglog( rbin, (rbin*parsec)*vphi_magbin)
	pl.xlim(3e-3, 3e0)
	pl.xticks(fontsize=24)
	pl.yticks(fontsize=24)
	#pl.minorticks_on()
	#pl.tick_params('both',length=8, width=1, which='minor')
	#pl.tick_params('both',length=10, width=1.5, which='major')
	pl.gcf().subplots_adjust(bottom=0.15)
	pl.ylabel('Specific Ang momentum', fontsize=25)
	pl.xlabel('$r$ $({\\rm pc})$', fontsize=25)
	pl.rc('text', usetex=False)
	#pl.title('Particle {0} {1:03d} kyr after formation'.format(particle_name, int(partAge/1e3)))

def total_mass_plot(rbin, mTbin):
	pl.clf()
	pl.rc('text', usetex=True)
	pl.rc('font', family='serif')
	pl.loglog( rbin, mTbin)
	#pl.axhline(chosen_ratio_number*particleMass, color = 'g')
#	if not (withNoParticle) or (withdisk):
#		if disk_radius_index != 0 or particleMass != 0.0:
#			x_plot = [r_star, r_star]
#			y_plot = [1e32, chosen_ratio_number * disk_particle_mass]
#			pl.plot (x_plot, y_plot, '.--', color = 'g')
#			pl.plot (r_star, chosen_ratio_number * disk_particle_mass, color = 'g', marker = 'o')
#			pl.annotate('$R_*$ = ' + "%.2f"%r_star + 'pc', xy=(r_star, 1e32), xytext=(r_star, 1.25e32))
#
	pl.ylim(1e32, 1e37)
	pl.xlim(3e-3, 3e0)
	pl.xticks(fontsize=24)
	pl.yticks(fontsize=24)
	#pl.minorticks_on()
	#pl.tick_params('both',length=8, width=1, which='minor')
	#pl.tick_params('both',length=10, width=1.5, which='major')
	pl.gcf().subplots_adjust(bottom=0.15)
	pl.ylabel(r'$M$ $({\rm g})$', fontsize=25)
	pl.xlabel(r'$r$ $({\rm pc})$', fontsize=25)
	pl.rc('text', usetex=False)
	#pl.title('Particle {0} {1:03d} kyr after formation'.format(particle_name, int(partAge)))
	#pl.title('Total Mass')
	
def Mdot_plot(rbin, vrbin, rhobin, mdotbin):
	pl.clf()
	pl.rc('text', usetex=True)
	pl.rc('font', family='serif')
	pl.loglog( rbin, -4.0*pi*rhobin*vrbin*(parsec*rbin)**2*sec_in_yr/Msun) 
#	if not (withNoParticle) or (withdisk):
#		if disk_radius_index != 0 or particleMass != 0.0:
#			x_plot = [r_star, r_star]
#			y_plot = [1e-6, 1e-2]
#			pl.plot (x_plot, y_plot, '.--', color = 'g')
#			pl.plot (r_star, chosen_ratio_number * disk_particle_mass, color = 'g', marker = 'o')
#			pl.annotate('$R_*$ = ' + "%.2f"%r_star + 'pc', xy=(r_star, 1e32), xytext=(r_star, 1.25e32))
#
	pl.xlim(3e-3, 3e0)
	pl.ylim(1e-6, 1e-2)
	pl.xticks(fontsize=24)
	pl.yticks(fontsize=24)
	#pl.minorticks_on()
	#pl.tick_params('both',length=8, width=1, which='minor')
	#pl.tick_params('both',length=10, width=1.5, which='major')
	pl.gcf().subplots_adjust(bottom=0.15)
	pl.ylabel(r'$\dot{M}$ $({\rm M_{\odot}\, yr}^{-1})$', fontsize=25)
	pl.xlabel(r'$r$ $({\rm pc})$', fontsize=28)
	pl.rc('text', usetex=False)
	#pl.title('Particle {0} {1:03d} kyr after formation'.format(particle_name, int(partAge)))
	#pl.title('$\dot{M}$')

def velocity_vs_density_plot(rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, vphi_magbin, rhobin):
	pl.clf()
	pl.rc('text', usetex=True)
	pl.rc('font', family='serif')

	#pl.loglog( rhobin/mp, -vrbin/1e5, 'bv-', label='Infall Velocity')
	pl.loglog( rhobin/mp, vrmsbin/1e5, 'g.-', label= '$v_T$')#  r'\textbf{time} (s)'
	#pl.loglog( rhobin/mp, vKbin/1e5, 'r--', label='Keplerian Velocity')  
	pl.loglog( rhobin/mp, vphi_magbin/1e5, 'k+', label='$v_\phi$')
	murray_chang_velocity = np.sqrt(vphi_magbin*vphi_magbin + vrmsbin*vrmsbin)
	pl.loglog( rhobin/mp, murray_chang_velocity/1e5,'.-.', color='orange', label='$\sqrt{v_T^2 + v_\phi^2}$')
	#pl.loglog( rhobin/mp, sum_of_velocities/1e5, 'c', label='Sum Velocity')
	pl.annotate('M* = ' + "%.2f"% float(particleMass/Msun) + 'Msun', xy=(1e7, 1e-1), xytext=(1e6, 1e-1))
	pl.legend(loc=0, fontsize='medium', frameon=False)
	pl.xlim(1e4, 1e7)
	pl.ylim(1e-1, 5e0)
	pl.xticks(fontsize=24)
	pl.yticks(fontsize=24)
	#pl.minorticks_on()
	#pl.tick_params('both',length=8, width=1, which='minor')
	#pl.tick_params('both',length=10, width=1.5, which='major')
	pl.gcf().subplots_adjust(bottom=0.15)
	pl.ylabel(r'$v_T, \, v_\phi, \, \sqrt{v_T^2 + v_\phi^2}$ $({\rm km\, s}^{-1})$', fontsize=25)
	pl.xlabel('$N$ $({\\rm cm}^{-3})$', fontsize=25)
	pl.rc('text', usetex=False)
	#pl.title('Particle {0} {1:03d} kyr after formation'.format(particle_name, int(partAge)))
	#pl.title('Velocity by Density')

def Toomre_Q(rbin, rhobin, vrbin, vrmsbin, vphi_magbin):
	Q = np.zeros(len(rbin))
	Q = (vphi_magbin * vrmsbin) / (4.* np.pi * G * (rbin*parsec)**2 * rhobin)
	#print Q
	pl.clf()
	pl.rc('text', usetex=True)
	pl.rc('font', family='serif')

	pl.plot( rbin, Q, '.-', label='$Q$')  
	pl.axhline(1.0, color = 'black')  
	#pl.annotate('M* = ' + "%.2f"% float(particleMass/Msun) + 'Msun', xy=(1e7, 1e-1), xytext=(1e6, 1e-1))
	#pl.legend(loc=0, fontsize='medium', frameon=False)
	pl.xscale('log')
	pl.xlim(3e-3, 3e-2)
	pl.ylim(1e-1, 5e0)

	pl.xticks(fontsize=24)
	pl.yticks(fontsize=24)
	#pl.tick_params('both',length=8, width=1, which='minor')
	#pl.tick_params('both',length=10, width=1.5, which='major')
	pl.gcf().subplots_adjust(bottom=0.15)
	pl.ylabel('$Q$', fontsize=25)
	pl.xlabel('$r$ $({\\rm pc})$', fontsize=25)
	pl.rc('text', usetex=False)
	#pl.title('Particle {0} {1:03d} kyr after formation'.format(particle_name, int(partAge)))


def Pressure_Gravity_ratio(rbin, rhobin, vrmsbin, mTbin):
	global Big_sum
	global Call_count
	dr =  np.zeros(len(rbin))
	Gravity_bin =  np.zeros(len(rbin))
	P_thermal_bin =  np.zeros(len(rbin))
	P_turbulent_bin =  np.zeros(len(rbin))
	P_rotational_bin =  np.zeros(len(rbin))

	dP_thermal =  np.zeros(len(rbin))
	dP_turbulent =  np.zeros(len(rbin))
	dP_rotational =  np.zeros(len(rbin))
	Thermal_ratio =  np.zeros(len(rbin))
	Turbulent_ratio =  np.zeros(len(rbin))
	Rotational_ratio =  np.zeros(len(rbin))

	#Set rbin to be in cm
	rbin_cgs = rbin * parsec
	
	P_thermal_bin = rhobin * c_s*c_s
	P_turbulent_bin = rhobin * vrmsbin*vrmsbin
	P_rotational_bin = rhobin * vphi_magbin*vphi_magbin

	dP_thermal[1:-1] = P_thermal_bin[2:] - P_thermal_bin[:-2]
	dP_turbulent[1:-1] = P_turbulent_bin[2:] - P_turbulent_bin[:-2]
	dP_rotational[1:-1] = P_rotational_bin[2:] - P_rotational_bin[:-2]
	dr[1:-1] = rbin_cgs[2:] - rbin_cgs[:-2]
	Gravity_bin = G * mTbin / (rbin_cgs*rbin_cgs)
	
	Thermal_numerator = dP_thermal / (dr*rhobin)
	Turbulent_numerator = dP_turbulent / (dr*rhobin)
	rotational_numerator = dP_rotational / (dr*rhobin)

	Thermal_ratio = -Thermal_numerator / Gravity_bin
	Turbulent_ratio = -Turbulent_numerator / Gravity_bin
	Rotational_ratio = -rotational_numerator / Gravity_bin
	Sum_of_ratios = Thermal_ratio + Turbulent_ratio

	Big_sum = Big_sum + Sum_of_ratios
	print Big_sum
	Call_count = Call_count + 1

	pl.clf()

	#pl.loglog(rbin, Thermal_ratio, color = 'black')
	#pl.loglog(rbin, Turbulent_ratio, color = 'r')
	pl.loglog(rbin, Rotational_ratio, '--', color = 'blue', label = '$P_{v_\phi}$')
	#pl.loglog(rbin, Big_sum/Call_count, color = 'r', label = Call_count)
	pl.loglog(rbin, Sum_of_ratios, '-', color = 'g', label = '$P_{c_s}$ + $P_T$')
	pl.legend(loc=0, fontsize='22', frameon=False)
	pl.axhline(1.0, color='black')
	pl.xlim(3e-3, 3e0)
	pl.ylim(1e-2, 1e1)
	pl.xticks(fontsize=24)
	pl.yticks(fontsize=24)
	#pl.minorticks_on()
	#pl.tick_params('both',length=8, width=1, which='minor')
	#pl.tick_params('both',length=10, width=1.5, which='major')
	#pl.gcf().subplots_adjust(bottom=0.15)
#	pl.rc('text', usetex=True)
#	pl.rc('font', family='serif')
	pl.ylabel(r'${(dP/dr)}/{\rho g}, \, {v_{\phi}^2}/{rg}$', fontsize=25)
	pl.xlabel('$r$ $({\\rm pc})$', fontsize=25)
	#pl.rc('text', usetex=False)
	#pl.title('Ratio of Thermal+Turbulent Pressure Gradient to Gravity')
	#pl.title('Particle {0} {1:03d} kyr after formation'.format(particle_name, int(partAge)))

def find_r_star(rbin, mTbin, particleMass, r_star):
##Find r_*. Since r is ordered, I can use bisect to locate the index 
#	where m_t= 2 m_*
#	index_r_star = bisect.bisect_left(mass, 1.95 * mass[0])
#	r_star = r[index_r_star]
#
	global r_star1
	print "Now finding R * "
	current_closest1 = 0.05
	current_closest2 = 0.95
	first_one = 0

	wanted_ratio = 1 / chosen_ratio_number
	lower_end = wanted_ratio - .10
	upper_end = wanted_ratio + .10
	x1 = False
	x0 = False

	for i in range(len(rbin)):
		if mTbin[i] != 0.0:
			mass_ratio = disk_particle_mass / mTbin[i]
		else :
			mass_ratio = 1.0
			continue
		#print mass_ratio
		#if mass_ratio > 0.45 and mass_ratio < 0.7:
		if mass_ratio > wanted_ratio and mass_ratio < upper_end:
			#print mass_ratio, disk_particle_mass, mTbin[i]
			#if mass_ratio < current_closest2 and first_one == 0:
			if mass_ratio < current_closest2:
				current_closest2 = mass_ratio
				x1 = rbin[i]
				y1 = mTbin[i]

		#if mass_ratio > 0.2 and mass_ratio < 0.45:
		if mass_ratio > lower_end and mass_ratio < wanted_ratio:
			#print mass_ratio, disk_particle_mass, mTbin[i]
			#if mass_ratio > current_closest1:
			if mass_ratio > current_closest1 and first_one == 0:
				current_closest1 = mass_ratio
				x0 = rbin[i]
				y0 = mTbin[i]
				first_one = 1

	if x1 == False:
		upper_end = wanted_ratio + .25
		print 'Upper bound was not satisfied'
		for i in range(len(rbin)):
			if mTbin[i] != 0.0:
				mass_ratio = disk_particle_mass / mTbin[i]
			else :
				mass_ratio = 1.0
				continue
			if mass_ratio > wanted_ratio and mass_ratio < upper_end:
				#print mass_ratio, disk_particle_mass, mTbin[i]
				if mass_ratio < current_closest2:
					current_closest2 = mass_ratio
					x1 = rbin[i]
					y1 = mTbin[i]

	if x0 == False:
		lower_end = wanted_ratio - .25
		first_one = 0
		print 'Lower bound was not satisfied'
		for i in range(len(rbin)):
			if mTbin[i] != 0.0:
				mass_ratio = disk_particle_mass / mTbin[i]
			else :
				mass_ratio = 1.0
				continue
			if mass_ratio > lower_end and mass_ratio < wanted_ratio:
				#print mass_ratio, disk_particle_mass, mTbin[i]
				if mass_ratio > current_closest1 and first_one == 0:
					current_closest1 = mass_ratio
					x0 = rbin[i]
					y0 = mTbin[i]

	# We now have x0, y0 and x1, y1, time to interpolate and then find r_*
	#y = 4.0 * disk_particle_mass
	y = chosen_ratio_number * disk_particle_mass
	print y / disk_particle_mass
	r_star1 = x0 + (x1 - x0) * (y - y0) / (y1 - y0)
	print r_star1
	#print 2.0 * disk_particle_mass
	#return r_star1

def r_star_by_gasmass(rbin, mTbin, particleMass, r_star):
	global r_star1
	#Find the sphere of influence of the star by setting
	# the stellar + disk mass = gas mass - disk mass
	star_gravity = G * (particleMass + diskMass)
	gas_gravity = G * (mTbin - diskMass)


def density_derivative_profile(rbin, rhobin):
	# Find disk radius by derivative of density
	global dlnrho_dlnr
	global disk_radius_index
	
	outside_radius_set = False
	disk_radius_index = 0
	drho = np.zeros(len(rbin))
	dr = np.zeros(len(rbin))
	dlnrho_dlnr = np.zeros(len(rbin))
	ratio_deriv = np.zeros(len(rbin))

	drho[1:-1] = rhobin[2:] - rhobin[:-2] 
	dr[1:-1] = rbin[2:] - rbin[:-2] 
	ratio_deriv[1:-1] = drho[1:-1]/dr[1:-1]
	ratio_deriv = ratio_deriv * -1.0
	dlnrho_dlnr = ratio_deriv * rbin / rhobin
	for i in range(len(rbin), -1, -1):
		if i != 0:
			if i != len(rbin):
				try:
					math.log10(dlnrho_dlnr[i])
				except ValueError:
					print 'Attempted to log a negative number at bin ', i
					continue
				#print dlnrho_dlnr[i], rbin[i]
				if outside_radius_set == False:
					if dlnrho_dlnr[i] <=0.99 and rbin[i] <= 0.1:
						outside_radius_set = True
						disk_radius_index = i
						print 'current disk_radius is', disk_radius_index
						print 'corrsponds to rbin of: ', rbin[i]
				

	if outside_radius_set == False:
		print 'No disk'

def disk_radius_by_vphi(rbin, vrbin, vrmsbin, vphi_magbin):
	global disk_radius_index
	boolean_r_mask = vphi_magbin > 0.99*vrbin
	boolean_t_mask = vphi_magbin > 0.99*vrmsbin
	if np.where(boolean_r_mask)[0].size:
		last_r_index = np.where(boolean_r_mask)[0][-1]
	elif not np.where(boolean_r_mask)[0].size:
		last_r_index = 69
	if np.where(boolean_t_mask)[0].size:
		last_t_index = np.where(boolean_t_mask)[0][-1]
	elif not np.where(boolean_t_mask)[0].size:
		last_t_index = 69
	print last_r_index
	print last_t_index
	if last_r_index < last_t_index:
		disk_radius_index = last_r_index
	elif last_t_index < last_r_index:
		disk_radius_index = last_t_index
	elif last_t_index == last_r_index and last_r_index != 69:
		disk_radius_index = last_r_index
	elif last_t_index == last_r_index and last_r_index == 69:
		disk_radius_index = 0
		print 'Appears there is NO rotationally supported disk'
	print 'This is rbin[disk_radius_index] in pcs', rbin[disk_radius_index]
#	radius_disk_bin = rbin
#	print radius_disk_bin[boolean_t_mask]
#	print radius_disk_bin[boolean_r_mask]
#	print radius_disk_bin[boolean_r_mask][boolean_t_mask]
#	radius_disk_bin = radius_disk_bin[boolean_r_mask][boolean_t_mask]
#	r_disk = radius_disk_bin[-1]
#	if r_disk == rbin[disk_radius_index]:
#		print 'This is the disk radius in pc: ', r_disk
#	elif r_disk != rbin[disk_radius_index]:
#		print 'We"ve run into trouble with finding the disk radius'
#		sys.exit()

def naming_convention(quad_number, particle_number):
	#global particle_name
	particle_name = 'Q' + str(quad_number) + '_' + str(particle_number)
	print particle_name
	return particle_name

def pdf_close(pdf_timestep):
	pdf_timestep.close()

def png_save(fileout):
	print 'saving file: ', fileout
	pl.savefig(fileout)

def disperse_2_plotters(quad, framestep, compare_files, particle_number, output_format, rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, rhobin, mdotbin, norm, angXbin, angYbin, angZbin, vphi_magbin, sum_of_velocities, partMass, part_creation_time, current_time, vrms_r_bin, vrms_l_bin, vrms_theta_bin, vrms_phi_bin, all_independent, velocity_alone, Toomre_Q_alone, density_alone, angmv_alone, mass_alone, mdot_alone, veldens_alone, Pressure_alone):
	
	if (all_independent):
		velocity_alone = True
		if not (withNoParticle) or (withdisk):
			Toomre_Q_alone = True
		density_alone = True
		angmv_alone = True
		mass_alone = True
		mdot_alone = True
		veldens_alone = True
		Pressure_alone = True
	if (velocity_alone):
		fileout = "velocity_{0}_{1:04d}_{2}{3}.{4}".format(quad, framestep, compare_files, particle_number, output_format)
		velocity_plot(rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, vphi_magbin)
		png_save(fileout)
		fileout = "Turbulent_{0}_{1:04d}_{2}{3}.{4}".format(quad, framestep, compare_files, particle_number, output_format)
		Turbulent_velocity_plot(rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, vphi_magbin)
		png_save(fileout)
	if (Toomre_Q_alone):
		fileout = "Toomre_Q_{0}_{1:04d}_{2}{3}.{4}".format(quad, framestep, compare_files, particle_number, output_format)
		Toomre_Q(rbin, rhobin, vrbin, vrmsbin, vphi_magbin)				
		png_save(fileout)
	if (veldens_alone):
		fileout = "vel_dens_{0}_{1:04d}_{2}{3}.{4}".format(quad, framestep, compare_files, particle_number, output_format)
		velocity_vs_density_plot(rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, vphi_magbin, rhobin)
		png_save(fileout)
	if (density_alone):
		fileout = "density_{0}_{1:04d}_{2}{3}.{4}".format(quad, framestep, compare_files, particle_number, output_format)
		density_plot(rbin, mTbin, rhobin)
		png_save(fileout)
	if (Pressure_alone):
		fileout = "Pressure_{0}_{1:04d}_{2}{3}.{4}".format(quad, framestep, compare_files, particle_number, output_format)
		Pressure_Gravity_ratio(rbin, rhobin, vrmsbin, mTbin)
		png_save(fileout)
	if (angmv_alone):
		fileout = "angmomentum_{0}_{1:04d}_{2}{3}.{4}".format(quad, framestep, compare_files, particle_number, output_format)
		ang_moment_plot(rbin, mTbin, vphi_magbin)
		png_save(fileout)
	if (mass_alone):
		fileout = "mass_{0}_{1:04d}_{2}{3}.{4}".format(quad, framestep, compare_files, particle_number, output_format)
		total_mass_plot(rbin, mTbin)
		png_save(fileout)
	if (mdot_alone):
		fileout = "mdot_{0}_{1:04d}_{2}{3}.{4}".format(quad, framestep, compare_files, particle_number, output_format)
		Mdot_plot(rbin, vrbin, rhobin, mdotbin)
		png_save(fileout)

def least_squares_fitting(log_x_sorted, log_y_sorted):
#
	#A = n.vstack([l_Mass_total_sort, n.ones(len(l_Mass_total_sort))]).T
	#m, c = n.linalg.lstsq(A, l_L_to_M_sort)[0]
	#L_to_M_fit = 10.**(m*l_Mass_total_sort + c)
	A = np.vstack([log_x_sorted, np.ones(len(log_x_sorted))]).T
	m, c = np.linalg.lstsq(A, log_y_sorted)[0]
	x_to_y_fit = 10.**(m*log_x_sorted + c)
	print('slope =', m)
	return x_to_y_fit
#

def getValues( particleFile, value) :
        names = np.char.strip(particleFile["particle names"].value)
        starParticles = particleFile["tracer particles"].value
        ind = (np.char.strip(names) == value)[:,0]
        return starParticles[:, ind].flatten()

def File_opener(filename):
    particleFile = h5py.File(filename)
#    print particleFile.keys()
    # Get all the star particles first
    if("tracer particles" in particleFile.keys()) :
	    return particleFile
    else :
	    particleFile.close()
	    return None

def obtain_Particle_ID(filename) :
	#filename = "quad{0}/{1}{2:04d}".format(quad,part_prefix, fileIndex) 
	#filename = '{0}_{1}_{2:04d}'.format(part_prefix, quad, framestep))

	file_exist = glob.glob(filename)
	if not file_exist:
		print 'Need the particle file to plot. Not just the output file.'
	particleFile = File_opener(filename)
	if( particleFile == None) : 
		partIndices = [] #an empty list 0
	if( particleFile != None) : 
		#current_time = particleFile['real scalars'][0][1]
		partIndices = getValues(particleFile, "tag").astype("int")
		#print partIndices
		#print type(partIndices)
		#sys.exit()
		#partCreateTime = getValues(particleFile, "creation_time")
		particleFile.close()
#	print partIndices
	print len(partIndices)
	return partIndices


parser = argparse.ArgumentParser(description = " Requires: start number, end number, stride_length, file_name, which_plots")

parser.add_argument('start', metavar='N1', type=int, help='File start Point')
parser.add_argument('end', metavar='N2', type=int, help='File end Point (does not include this file)')
parser.add_argument('step', metavar='N3', type=int, help='Step size')
parser.add_argument('ParticleID', metavar='N4', type=int, nargs='?', default=42, help='Particle ID you want to reduce.')
parser.add_argument('--bigsphere', action='store_true')
parser.add_argument('--smallsphere', action='store_true')
parser.add_argument('--particle', action='store_true')
parser.add_argument('--noparticle', action='store_true')
parser.add_argument('--shell', action='store_true')
parser.add_argument('--shellsphere', action='store_true')
parser.add_argument('--velocity', action='store_true')
parser.add_argument('--density', action='store_true')
parser.add_argument('--veldens', action='store_true')
parser.add_argument('--mdot', action='store_true')
parser.add_argument('--mass', action='store_true')
parser.add_argument('--angmv', action='store_true')
parser.add_argument('--allparticles', action='store_true')
parser.add_argument('--allplots', action='store_true')
parser.add_argument('--disk', action='store_true')
parser.add_argument('--Q', action='store_true')
parser.add_argument('--pressure', action='store_true')
parser.add_argument('--pdf', action='store_true')

args = parser.parse_args()

withBigSphere = args.bigsphere
withSmallSphere = args.smallsphere
withParticle = args.particle
withNoParticle = args.noparticle
withShell = args.shell
withShellSphere = args.shellsphere
withdisk = args.disk
withParticleIDValue = args.ParticleID
velocity_alone = args.velocity
veldens_alone = args.veldens
density_alone = args.density
angmv_alone = args.angmv
mdot_alone = args.mdot
mass_alone = args.mass
Pressure_alone = args.pressure
Toomre_Q_alone = args.Q
all_independent = args.allplots
withAllParticles = args.allparticles
withPDF = args.pdf

mp = 1.6e-24
pi = np.pi
parsec = 3e18
sec_in_yr = np.pi* 1e7
Msun = 2e33
G = 6.67e-8
c_s = 2.64e4

Gravity_Turned_On = 6.317465983e14

# These have to do with calculating the sphere of influence
r_star = 0.0
chosen_ratio_number = 3.0
particleMass = 0.0

#Attempt at a hack of an avg
Call_count = 0

Big_sum = np.zeros(70)
if (withBigSphere):
	compare_files = 'bigsphere_'
	vel_Plot_label = '3.0 pc'
elif (withSmallSphere):
	compare_files = 'smallsphere_'
	vel_Plot_label = '0.01 pc'
elif (withParticle):
	compare_files = 'part_'
	vel_Plot_label = 'Bulk by Particle Vel'	
elif (withShell):
	compare_files = 'shell_'
	vel_Plot_label = 'Bulk by Shell'
elif (withShellSphere):
	compare_files = 'shellsphere_'
	vel_Plot_label = 'Bulk velocity removed by Sphere in Shell'

#if (withNoParticle):
#	vel_Plot_label = 'No Particles'	
#
if (withPDF):
	output_format = 'pdf'
else:
	output_format = 'png'

file_prefix = 'rad_profile'
part_prefix = 'BB_hdf5_part'
quad = os.getcwd()[-5:]
quad_number = quad[-1]

# The file looks like this:
#'{file_prefix}{framestep}_{compare_files}_{particle_number}.out'
# i.e. rad_profile_0218_part_0000.out

for framestep in range(args.start,args.end,args.step) :   
	# For each file in the given set, check to see if the reduced output for it exists.
	file_exist = glob.glob('{0}_{1}_{2:04d}_{3}*.out'.format(file_prefix, quad, framestep, compare_files))
	if not file_exist:
		print 'File: "{0}_{1}_{2:04d}_{3}*.out" does not exist!'.format(file_prefix, quad, framestep, compare_files)
		file_exist = glob.glob('{0}_{1}_{2:04d}_{3}*.out'.format(file_prefix, quad, framestep, compare_files))
		if not file_exist:
			print 'File: "{0}_{1}_{2:04d}_{3}*.out" does not exist!'.format(file_prefix, quad, framestep, compare_files)
			continue
	# If the file exists, 
	part_filename = '{0}_{1:04d}'.format(part_prefix, framestep)
	print part_filename
	#particle_list, current_time = obtain_Particle_ID(filename)
	particle_list = obtain_Particle_ID(part_filename)
#	print len(particle_list)
	if len(particle_list) == 0:
		print 'Particles have not been created yet'
		print 'Setting the particle list to include the specified particle id.'
		particle_list = [withParticleIDValue]
	print 'The total number of star particles in this timestep is: ', len(particle_list)
	print particle_list
	for j in range(len(particle_list)):
		particle_number = particle_list[j]
		if int(particle_number) == withParticleIDValue or (withAllParticles):
			print "Plotting particle", j+1, 'of', len(particle_list)
#			filein = '{0}_{1}_{2:04d}_{3}{4}.out'.format(file_prefix, quad, framestep, compare_files, particle_number)
			filein = '{0}_{1}_{2:04d}_{3}{4}.0.out'.format(file_prefix, quad, framestep, compare_files, particle_number)
			particle_file_exist =  glob.glob(filein)
			print filein
			if not particle_file_exist:
				print 'File: ', filein, ' does not exist!'
				continue
			# If the file in question exists, we then unpack
			rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, rhobin, mdotbin, norm, angXbin, angYbin, angZbin, vphi_magbin, sum_of_velocities, partMass, part_creation_time, current_time, vrms_r_bin, vrms_l_bin, vrms_theta_bin, vrms_phi_bin = np.loadtxt(filein, unpack=True)

			part_creation_time = part_creation_time[0]
			print part_creation_time
			current_time = current_time[0] # This can also be obtained from part_file.
			print 'current time is: ', current_time
			print 'time since gravity (kyr): ', (current_time - Gravity_Turned_On) / (3.14e10)
			partAge = current_time - part_creation_time
			partAge = partAge / (sec_in_yr * 1e3)
			if part_creation_time == 0.0: # withNoParticle:
				partAge = 0.0
			print 'particle age in thousands of years', partAge
			# Get the name of the particle for the title of the plots.
			#particle_name = naming_convention(quad_number, particle_number)
			#print partMass
			# convert to cgs, mTbin is already in cgs
			partMass = partMass * Msun

			# If only one particle set =, if a list, pick which particle its on
			if isinstance(partMass, float) ==True:
				particleMass = partMass
			else:
				particleMass = partMass[0]

			# Calculate the derivative of the density
			# And calculate the disk radius
			# Old method of finding disk radius
			#density_derivative_profile(rbin, rhobin)

			if not (withNoParticle) or (withdisk):
				disk_radius_by_vphi(rbin, vrbin, vrmsbin, vphi_magbin)
				diskMass = mTbin[disk_radius_index] - particleMass
			else:
				disk_radius_index = 0.0
				diskMass = 0.0
			# mTbin[i] is all the mass in the sphere inside, not just a shell

			print 'Disk Mass is: ', diskMass, 'g or ', diskMass/Msun, 'M_sun'
			print 'Star particle Mass is: ', particleMass, 'g or ', particleMass/Msun, '$M_sun$'
			try:
				print 'Disk/Star Mass Ratio is: ', diskMass / particleMass
			except ValueError:
				print 'No stellar mass so dividing by zero'
			disk_particle_mass = particleMass + diskMass

			if not (withNoParticle) or (withdisk):
				if disk_radius_index != 0 or particleMass != 0.0:
					try:
						find_r_star(rbin, mTbin, particleMass, r_star)
					except UnboundLocalError:
						r_star1 = 0.0
					# Now set r_star = to the global r_star1
					r_star = r_star1
					print 'R_* is: ', r_star, 'pcs or ', r_star * parsec, 'cms'
			try:
				disperse_2_plotters(quad, framestep, compare_files, particle_number, output_format, rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, rhobin, mdotbin, norm, angXbin, angYbin, angZbin, vphi_magbin, sum_of_velocities, partMass, part_creation_time, current_time, vrms_r_bin, vrms_l_bin, vrms_theta_bin, vrms_phi_bin, all_independent, velocity_alone, Toomre_Q_alone, density_alone, angmv_alone, mass_alone, mdot_alone, veldens_alone, Pressure_alone)
			except ValueError:
				print "Can not plot this particle"
		else:
			#print "The particle ID's don't match, proceeding to next particle."
			continue
