"""
This script is similar to FLASH_plot.py eventual intention is to combine the two into one script.
This script obtains its data from the rad_profile*.out produced by all_profiles*.py as well as the RAMSES sink file for the star particle data.
"""
######## All import statements #################
import matplotlib.pyplot as plt
import numpy as np
import argparse
import glob
import sys
import os
import math
import warnings # To deal with empty sinkfiles.
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import matplotlib
from matplotlib import rcParams
matplotlib.use("Agg")
warnings.filterwarnings('error','loadtxt: Empty input file:')

def Setup_plot_window():
	global Init_matplotParams
	if Init_matplotParams == False:
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
		#xtick.minor.width    : 0.5    # minor tick width in points
		Init_matplotParams = True

	plt.clf()
	if (args.rstar):
		if (disk_radius_index != 0) or (particleMass != 0.0) :
			x_plot = [r_star, r_star]
			y_plot = [1e-2, 1e4]
			plt.plot (x_plot, y_plot, '.--', color = 'g')
			plt.annotate('$R_*$ = ' + "%.2f"%r_star + 'pc', xy=(r_star, 1e1), xytext=(r_star, 1.25e-2), fontsize=20)
	#plt.legend(loc=3, fontsize=21, frameon=False, ncol=2)
	#plt.title('{0} framestep {1:04d}'.format(vel_Plot_label, framestep))
	#plt.title('Particle {0} {1:3d} kyr after formation'.format(particle_name, int(partAge)))
	plt.xticks(fontsize=24)
	plt.yticks(fontsize=24)
	plt.minorticks_on()
	plt.tick_params('both',length=8, width=1, which='minor')
	plt.tick_params('both',length=10, width=1.5, which='major')
#	plt.gcf().subplots_adjust(bottom=0.15)
	plt.tight_layout()


def Velocity_plot(rbin, vrbin, vrmsbin, vmagbin, vKbin, vphi_magbin, sum_velocities, MC15_velocity):
#	if( args.avg) :
#		host = host_subplot(111, axes_class=AA.Axes)
#		host.set_xscale("log", nonposx='clip')
#		host.set_yscale("log", nonposy='clip')
#		print rbin, vrbin/1e5
#		#sys.exit()
#		host.errorbar(rbin, -vrbin/1e5, xerr=0.0*rbin, yerr=sigma, color='bv-', label='$u_r$')
#		host.loglog( rbin, vrmsbin/1e5, 'g.-', label='$v_T$')  
#		host.loglog( rbin, vKbin/1e5, 'r--', label='$v_K$')  
#		host.loglog( rbin, vphi_magbin/1e5, 'k+', label='$v_\phi$')
#		host.axhline( c_s/1e5, color = 'black',label='$c_s$')  
#		#host.ylim(5e-3, 3e1)
#		#host.xlim(3e-3, 3e0)
#		host.ylabel(r'$u_r, \, v_T, \, v_K, \, v_\phi, \, c_s$ $({\rm km\, s}^{-1})$', fontsize = 25)
#		host.xlabel(r'$r$ $({\rm pc})$', fontsize=25)
	plt.loglog( rbin, -vrbin/1e5, 'bv-', label='$u_r$')
	plt.loglog( rbin, vrmsbin/1e5, 'g.-', label='$v_T$')  
	#plt.loglog( rbin, vmagbin/1e5, 'm.', label='$v_{Total}$')  
	plt.loglog( rbin, vKbin/1e5, 'r--', label='$v_K$')  
	#plt.loglog( rbin, MC15_velocity/1e5, '.--', color='orange', label='v_t + v_phi')  
	plt.loglog( rbin, vphi_magbin/1e5, 'k+', label='$v_\phi$')
	#plt.loglog( rbin, sum_velocities/1e5, 'c', label='$\sqrt{v_r^2 + v_T^2 + v_\phi^2}$')  
	# Need to mask this so that if not fully refined still get values.
	#power_law_fit = least_squares_fitting(log_x_sorted=rbin, log_y_sorted=vrmsbin)	
	#plt.loglog( rbin, power_law_fit, 'k.-', label='$powerlaw$')
	#print power_law_fit
	#sys.exit()
	plt.axhline( c_s/1e5, color = 'black',label='$c_s$')  
	plt.ylim(5e-3, 3e1)
	plt.xlim(3e-3, 3e0)
	plt.ylabel(r'$u_r, \, v_T, \, v_K, \, v_\phi, \, c_s$ $({\rm km\, s}^{-1})$', fontsize = 25)
	plt.xlabel(r'$r$ $({\rm pc})$', fontsize=25)
	
def Turbulent_velocity_plot(rbin, vrmsbin, vrms_r_bin, vrms_l_bin):
	plt.loglog( rbin, vrmsbin/1e5, 'go-', label='$v_T$', lw = 1)  
	plt.loglog( rbin, vrms_r_bin/1e5, '-', color='blue', label='$v_{T,r}$')  
	plt.loglog( rbin, vrms_l_bin/1e5, 'r--', label='$v_{T,l}$', lw = 2)  
	plt.axhline( c_s/1e5, color = 'black',label='$c_s$')  
	plt.legend(loc=0, fontsize=22, frameon=False)
	plt.ylim(1e-1, 3e0)
	plt.xlim(3e-3, 3e0)
	plt.ylabel(' $v_T, \\, v_{T,r}, \\, v_{T,l}$ $({\\rm km\\, s}^{-1})$', fontsize=25)
	plt.xlabel('$r$ $(\\rm pc)$', fontsize=25)
	
def Density_plot(rbin, rhobin):
	plt.clf()
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
	plt.xticks(fontsize=30)
	plt.yticks(fontsize=30)
	#par1.yticks(fontsize= 20)
	host.set_ylabel('$N$ $({\\rm cm^{-3}})$', fontsize=25)
	host.set_xlabel('$r$ $({\\rm pc})$', fontsize=25)
	par1.set_ylabel('$\\rho$ $({\\rm g \\, cm}^{-3})$', fontsize=25)
        host.axis["left"].label.set_fontsize(25)
        host.axis["bottom"].label.set_fontsize(25)
        par1.axis["right"].label.set_fontsize(25)
	# rhobin should be in g/cm**3
	host.loglog(rbin, rhobin/mp, label='Number density') 	
	#host.loglog(rbin, dlnrho_dlnr, label='derivative')
	#plt.legend()
	# If plotting multiple things, this may be useful
	#host.axis["left"].label.set_color(p1.get_color())
	#par1.axis["right"].label.set_color(p2.get_color())
	
def Ang_moment_plot(rbin, vphi_magbin):
	#plt.loglog( rbin, (rbin*parsec)*vphi_magbin/np.sqrt(G*mTbin*(parsec*rbin)))
	plt.loglog( rbin, (rbin*parsec)*vphi_magbin)
	plt.xlim(3e-3, 3e0)
	plt.ylabel('Specific Ang momentum', fontsize=25)
	plt.xlabel('$r$ $({\\rm pc})$', fontsize=25)

def Total_mass_plot(rbin, mTbin):
	plt.loglog( rbin, mTbin)
	#plt.axhline(chosen_ratio_number*particleMass, color = 'g')
	plt.ylim(1e32, 1e37)
	plt.xlim(3e-3, 3e0)
	plt.ylabel(r'$M$ $({\rm g})$', fontsize=25)
	plt.xlabel(r'$r$ $({\rm pc})$', fontsize=25)
	
def Mdot_plot(rbin, py_mdot, mdotbin):
#	plt.loglog( rbin, -4.0*pi*rhobin*vrbin*(parsec*rbin)**2*yr/Msun) 
	plt.loglog( rbin, py_mdot) 
	plt.loglog( rbin, mdotbin) 
	plt.xlim(3e-3, 3e0)
	plt.ylim(1e-6, 1e-2)
	plt.ylabel(r'$\dot{M}$ $({\rm M_{\odot}\, yr}^{-1})$', fontsize=25)
	plt.xlabel(r'$r$ $({\rm pc})$', fontsize=28)

def Velocity_vs_Density_plot(rhobin, vrbin, vrmsbin, vKbin, vphi_magbin, particleMass, MC15_velocity):
	plt.loglog( rhobin/mp, vrmsbin/1e5, 'g.-', label= '$v_T$')
	plt.loglog( rhobin/mp, vphi_magbin/1e5, 'k+', label='$v_\phi$')
	plt.loglog( rhobin/mp, MC15_velocity/1e5,'.-.', color='orange', label='$\sqrt{v_T^2 + v_\phi^2}$')
	#plt.loglog( rhobin/mp, -vrbin/1e5, 'bv-', label='$v_r$')
	#plt.loglog( rhobin/mp, vKbin/1e5, 'r--', label='v_{Kep}')  
	#plt.loglog( rhobin/mp, sum_of_velocities/1e5, 'c', label='Sum Velocity')
	plt.annotate('M* = ' + "%.2f"% float(particleMass/Msun) + 'Msun', xy=(1e7, 1e-1), xytext=(1e6, 1e-1))
	plt.xlim(1e4, 1e7)
	plt.ylim(1e-1, 5e0)
	plt.ylabel(r'$v_T, \, v_\phi, \, \sqrt{v_T^2 + v_\phi^2}$ $({\rm km\, s}^{-1})$', fontsize=25)
	plt.xlabel('$N$ $({\\rm cm}^{-3})$', fontsize=25)

def Toomre_Q(rbin, rhobin, vrmsbin, vphi_magbin):#rbin, rhobin, vrbin, vrmsbin, vphi_magbin):
	Q = np.zeros(len(rbin))
	Q = (vphi_magbin * vrmsbin) / (4.* np.pi * G * (rbin*parsec)**2 * rhobin)

	plt.plot( rbin, Q, '.-', label='$Q$')  
	plt.axhline(1.0, color = 'black')  
	#plt.annotate('M* = ' + "%.2f"% float(particleMass/Msun) + 'Msun', xy=(1e7, 1e-1), xytext=(1e6, 1e-1))
	plt.xscale('log')
	plt.xlim(3e-3, 3e-2)
	plt.ylim(1e-1, 5e0)
	plt.ylabel(r'$Q$', fontsize=25)
	plt.xlabel(r'$r$ $({\rm pc})$', fontsize=25)

def Pressure_Gravity_ratio(rbin, rhobin, vrmsbin, vphi_magbin, Thermal_P_ratio, Turb_P_ratio, Rot_P_ratio):
	Sum_of_ratios = Thermal_P_ratio + Turb_P_ratio
	#plt.loglog(rbin, Thermal_P_ratio, color = 'black')
	#plt.loglog(rbin, Turb_P_ratio, color = 'r')
	plt.loglog(rbin, Rot_P_ratio, '--', color = 'blue', label = '$P_{v_\phi}$')
	plt.loglog(rbin, Sum_of_ratios, '-', color = 'g', label = '$P_{c_s}$ + $P_T$')
	plt.legend(loc=0, fontsize='22', frameon=False)
	plt.axhline(1.0, color='black')
	plt.xlim(3e-3, 3e0)
	plt.ylim(1e-2, 1e1)
	plt.ylabel(r'${(dP/dr)}/{\rho g}, \, {v_{\phi}^2}/{rg}$', fontsize=25)
	plt.xlabel('$r$ $({\\rm pc})$', fontsize=25)

def png_save(fileout):
	print 'saving file: ', fileout
	plt.savefig(fileout)

#################################
#################################
#################################
#################################
def Pressure_calculation():
	bins = len(rbin)

	P_thermal =  np.zeros(bins)
	P_turbulent =  np.zeros(bins)
	P_rotational =  np.zeros(bins)
	Gravity =  np.zeros(bins)

	dr =  np.zeros(bins)
	dP_thermal =  np.zeros(bins)
	dP_turbulent =  np.zeros(bins)
	dP_rotational =  np.zeros(bins)

	Thermal_P_ratio =  np.zeros(bins)
	Turb_P_ratio =  np.zeros(bins)
	Rot_P_ratio =  np.zeros(bins)
	#Set rbin to be in cm
	rbin_cgs = rbin * parsec

	P_thermal = rhobin * c_s*c_s
	P_turbulent = rhobin * vrmsbin*vrmsbin
	P_rotational = rhobin * vphi_magbin*vphi_magbin
	Gravity = G * mTbin / (rbin_cgs*rbin_cgs)

	dr[1:-1] = rbin_cgs[2:] - rbin_cgs[:-2]
	dP_thermal[1:-1] = P_thermal[2:] - P_thermal[:-2]
	dP_turbulent[1:-1] = P_turbulent[2:] - P_turbulent[:-2]
	dP_rotational[1:-1] = P_rotational[2:] - P_rotational[:-2]

	Thermal_P_ratio = -dP_thermal / dr / rhobin / Gravity
	Turb_P_ratio = -dP_turbulent / dr / rhobin / Gravity
	Rot_P_ratio = -dP_rotational / dr / rhobin / Gravity
#	Sum_of_ratios = Thermal_P_ratio + Turb_P_ratio

	return Thermal_P_ratio, Turb_P_ratio, Rot_P_ratio

def Additional_Calculations():
	# use this to figure out what % of vkep vr infall is.
	vr_vk_ratio = abs(np.nan_to_num(vrbin) / np.nan_to_num(vKbin))
#	print vr_vk_ratio
#	print 'Minimum % vr is of vk: ', vr_vk_ratio.min()
	# velocity and Veldens use this.
	MC15_velocity = np.sqrt(vrmsbin*vrmsbin + vphi_magbin*vphi_magbin)
	py_mdot = -4.0*pi*rhobin*vrbin*(parsec*rbin)**2*yr/Msun
	return MC15_velocity, py_mdot

def Find_r_star(rbin, mTbin, particleMass):
##Find r_*. Since r is ordered, I can use bisect to locate the index 
#	where m_t= 2 m_*
#	index_r_star = bisect.bisect_left(mass, 1.95 * mass[0])
#	r_star = r[index_r_star]
#
	print "Now finding R * "
	wanted_ratio = 1 / gas_to_particle_ratio
	lower_bound = wanted_ratio - .05
	upper_bound = wanted_ratio + .05
	lower_mass_ratio = 0.05
	upper_mass_ratio = 0.95
	r1 = False
	r0 = False
	index = 0
	infinite_loop_check = 0
#	print len(rbin)
	while index <= len(rbin) - 1:
#		print 'mTbin',index, mTbin[index]
		if mTbin[index] != 0.0:
			mass_ratio = disk_particle_mass / mTbin[index]
		else :
			index = index + 1
			continue
		if mass_ratio > wanted_ratio and mass_ratio < upper_bound:
			if mass_ratio < upper_mass_ratio:
				# This is at smaller radii, recall ratio is particle mass/gas mass
				# We want the largest radii that this holds for.
				upper_mass_ratio = mass_ratio
				r0 = rbin[index]
				m0 = mTbin[index]
				print 'r0', r0, m0, mass_ratio
		elif mass_ratio > lower_bound and mass_ratio < wanted_ratio:
			if mass_ratio > lower_mass_ratio:
				lower_mass_ratio = mass_ratio
				r1 = rbin[index]
				m1 = mTbin[index]
				print 'r1', r1, m1, mass_ratio
		index = index + 1
		if index >= len(rbin) and ((r0 == False) or (r1 == False)):
			infinite_loop_check = infinite_loop_check + 1
			if infinite_loop_check > 10:
				print "There is no stellar sphere of influence."
				break
			if r1 == False:
				upper_bound = upper_bound + .05
				index = 0
				if infinite_loop_check == 1:
					print 'Upper bound was not satisfied'
			if r0 == False:
				lower_bound = lower_bound - .05
				index = 0
				if infinite_loop_check == 1:
					print 'Lower bound was not satisfied'
	# We now have r0, m0 and r1, m1, time to interpolate and then find r_*
	m = gas_to_particle_ratio * disk_particle_mass
	try:
		r_star = r0 + (r1 - r0) * (m - m0) / (m1 - m0)
	except UnboundLocalError:
		print 'One of the bounds was not set.'
		print 'There is no sphere of influence for this particle in this timestep.'
		r_star = 0.0
#	print 'r0', r0
#	print 'r1 minus r0', r1 - r0
#	print 'mass ratio ', (m - m0) / (m1 - m0)
#	print 'r_star', r_star
	return r_star
#
#def r_star_by_gasmass(rbin, mTbin, particleMass):
#	#Find the sphere of influence of the star by setting
#	# the stellar + disk mass = gas mass - disk mass
#	star_gravity = G * (particleMass + diskMass)
#	gas_gravity = G * (mTbin - diskMass)
#	return r_star
#

def Density_derivative_profile(rbin, rhobin):
	# Find disk radius by derivative of density
	global disk_radius_index
	disk_radius_index = 0

	outside_radius_set = False

	drho = np.zeros(len(rbin))
	dr = np.zeros(len(rbin))
	drhodr = np.zeros(len(rbin))
	dlnrho_dlnr = np.zeros(len(rbin))

	drho[1:-1] = rhobin[2:] - rhobin[:-2] 
	dr[1:-1] = rbin[2:] - rbin[:-2] 
	drhodr[1:-1] = -1.0 * drho[1:-1]/dr[1:-1]
	dlnrho_dlnr = drhodr * rbin / rhobin

	for i in range(len(rbin), -1, -1):
		if i == 0 or i == len(rbin):
			# Ghost cells
			continue
		try:
			math.log10(dlnrho_dlnr[i])
		except ValueError:
#			print 'Attempted to log a negative number at bin ', i
			continue
		#print dlnrho_dlnr[i], rbin[i]
		if outside_radius_set == False:
			if dlnrho_dlnr[i] <=0.99 and rbin[i] <= 0.1:
				outside_radius_set = True
				disk_radius_index = i
				print 'Disk radius by density is: ', rbin[disk_radius_index]
	if outside_radius_set == False:
		print 'No disk'

def Disk_radius_by_vphi(rbin, vrbin, vrmsbin, vphi_magbin):
	global disk_radius_index
	boolean_r_mask = vphi_magbin > 0.99*vrbin
	boolean_t_mask = vphi_magbin > 0.99*vrmsbin
	if np.where(boolean_r_mask)[0].size:
		last_r_index = np.where(boolean_r_mask)[0][-1]
	else :
		last_r_index = len(rbin) - 1
	if np.where(boolean_t_mask)[0].size:
		last_t_index = np.where(boolean_t_mask)[0][-1]
	else :
		last_t_index = len(rbin) - 1
	if last_r_index <= last_t_index:
		disk_radius_index = last_r_index
	else:
		disk_radius_index = last_t_index
	if disk_radius_index == len(rbin) -1: #69:
		print 'There is no disk by comparison of v_phi to v_r and v_T'
		disk_radius_index = 0
	print 'vphi > v_r inside',last_r_index
	print 'vphi > v_T inside',last_t_index

	print 'This is rbin[disk_radius_index] in pcs', rbin[disk_radius_index]

def Least_squares_fitting(log_x_sorted, log_y_sorted):

	#A = n.vstack([l_Mass_total_sort, n.ones(len(l_Mass_total_sort))]).T
	#m, c = n.linalg.lstsq(A, l_L_to_M_sort)[0]
	#L_to_M_fit = 10.**(m*l_Mass_total_sort + c)
	A = np.vstack([log_x_sorted, np.ones(len(log_x_sorted))]).T
	m, c = np.linalg.lstsq(A, log_y_sorted)[0]
	x_to_y_fit = 10.**(m*log_x_sorted + c)
	print('slope =', m)
	return x_to_y_fit

def standardize_File_number(input_File_number):
	output_File_number = "%05d"%input_File_number
	if len(output_File_number) >= 6:
		print 'This system assumes less than one million output files.'
		sys.exit()
	return output_File_number


def obtain_avg_value(tot_bin, particle_count, quantity_to_avg_bin):
	""" This function calculates the average and standard deviation.
	"""
	global Values_4_std
	global sigma
	# Keep track of the values in the avg to determine the standard deviation.
	Values_4_std.append(quantity_to_avg_bin)
	# Determine the new total.
	tot_bin = tot_bin + quantity_to_avg_bin
#	print tot_bin
	mean = tot_bin/particle_count
#	print 'mean', mean
	sigma = np.sqrt((Values_4_std - mean)**2/particle_count)
#	print len(Values_4_std), len(sigma), particle_count
#	print 'Values_4_std', Values_4_std
#	print 'sigma', sigma
	return tot_bin, mean


############################################
############################################
############################################
############################################
parser = argparse.ArgumentParser(description = " Requires: start number, end number, stride_length, file_name, which_plots")
parser.add_argument('start', metavar='N1', type=int, help ='Start value for hdf5 files')
parser.add_argument('end', metavar='N2', type=int, help='End Value for hdf5 files, note runs until End-1')
parser.add_argument('step', metavar='N3', type=int, help='Stride length')
parser.add_argument('ParticleID', metavar='N4', type=int, nargs='?', default=0, 
		    help='Particle ID you want to reduce, use 0 for no particle density peak.')
parser.add_argument('bulk_vel_method', metavar='N5', type=str, nargs='?', default='shellsphere',
		    help='method of removing the bulk motion of the gas. Options are: shellsphere, bigsphere, smallsphere, particle, shell.')
parser.add_argument('center_mass_4avg', metavar='N6', type=int, nargs='?', default=1,
		    help='Mass you want to average around, in solar masses. Default is 1.')

parser.add_argument('--noparticle', action='store_true')
parser.add_argument('--allparticles', action='store_true')
parser.add_argument('--avg', action='store_true')
#What do you want to plot?
parser.add_argument('--allplots', action='store_true')
parser.add_argument('--velocity', action='store_true')
parser.add_argument('--turbulent', action='store_true')
parser.add_argument('--density', action='store_true')
parser.add_argument('--veldens', action='store_true')
parser.add_argument('--mdot', action='store_true')
parser.add_argument('--mass', action='store_true')
parser.add_argument('--angmv', action='store_true')
parser.add_argument('--pressure', action='store_true')
parser.add_argument('--Q', action='store_true')
#Want to plot rstar?
parser.add_argument('--rstar', action='store_true')
parser.add_argument('--pdf', action='store_true')

args = parser.parse_args()
global Init_matplotParams
Init_matplotParams = False
global Values_4_std
Values_4_std = []
global sigma
sigma = []
file_prefix = 'rad_profile'

center_mass_4_avg = args.center_mass_4avg

mp = 1.6e-24
pi = np.pi
parsec = 3e18
yr = pi * 1e7
Msun = 2e33
G = 6.67e-8
c_s = 2.64e4

# This has to do with calculating the sphere of influence
gas_to_particle_ratio = 3.0
plotting_list = []
accepted_plotting_list = ['velocity', 'turbulent', 'toomreQ', 'veldens', 'density', 'angmomentum', 'pressure', 'mass', 'mdot']
if args.allplots:
	plotting_list = accepted_plotting_list
else:
	if ( args.velocity):
		plotting_list.append('velocity')
	if ( args.turbulent):
		plotting_list.append('turbulent')
	#if not (args.noparticle) or (args.rstar): #Not sure why I had this condition.
	if ( args.Q):
		plotting_list.append('toomreQ')
	if ( args.veldens):
		plotting_list.append('veldens')
	if ( args.density):
		plotting_list.append('density')
	if ( args.angmv):
		plotting_list.append('angmomentum')
	if ( args.pressure):
		plotting_list.append('pressure')
	if ( args.mass):
		plotting_list.append('mass')
	if ( args.mdot):
		plotting_list.append('mdot')
if plotting_list == []:
	print "You have not selected the plots to create."
	print "Possible plots include: ", accepted_plotting_list
	sys.exit()

bulk_vel_accepted_strings = {'shellsphere', 'bigsphere', 'smallsphere', 'particle', 'shell'}
if args.bulk_vel_method in bulk_vel_accepted_strings:
	compare_file = args.bulk_vel_method
else:
	print "You have not selected a the method used to remove the bulk velocity."
	print "Possible methods include: ", bulk_vel_accepted_strings
	sys.exit()

if (args.pdf):
	output_format = 'pdf'
else:
	output_format = 'png'

if (args.avg):
	tot_vr = 0
	tot_vrms = 0
	tot_density = 0
	tot_mdot = 0
	deltaMass = center_mass_4_avg * 0.25
	particle_count = 0
#	if deltaMass > 1.0:
#		deltaMass = 1.0 # range should be no more than a solar mass

# The file looks like this:
#'{file_prefix}{framestep}_{compare_file}_{particle_number}.out'
# i.e. rad_profile_0218_part_0000.out

for timestep in range(args.start,args.end,args.step) :   
	# For each timestep, check to see if the info and sink files exist.
	File_number = standardize_File_number(timestep)
	infofile = "info_{0}.txt".format(File_number)
	sinkfile = "sink_{0}.info".format(File_number)
	if( not os.path.isfile( sinkfile)):# or ( not os.path.isfile( infofile)):
		#print "Either the sink or info file for timestep", File_number, "is missing."
		print "The sink file for timestep", File_number, "is missing."
		continue #return
	try :
		particle_list = np.loadtxt( sinkfile, usecols=[0], unpack=True, skiprows=3, comments="=")
	except UserWarning :
		print 'Sink file is empty'
		if (args.noparticle):
			particle_list = np.array([int(args.ParticleID)])
		else :
			print 'The sinkfile is empty and you have not requested to plot max density location.'
			print "If you want to, use the flag --noparticle."
			continue #return
	if particle_list.size == 1:
		particle_list  = np.array([particle_list])
	print "Plotting Timestep:", File_number
#	for j in np.nditer(particle_list):
	for j in range(len(particle_list)):
		particle_ID = int(particle_list[j])
		if particle_ID == args.ParticleID or (args.allparticles): # (withAllParticles):
			# Does the particle have a reduced output for this timestep?
			filein = '{0}_{1}_{2}_{3}.out'.format(file_prefix, File_number, compare_file, particle_ID)
			particle_file_exist =  glob.glob(filein)
			if not particle_file_exist:
				filein2 = '{0}_{1}_{2}_{3}.0.out'.format(file_prefix, File_number, compare_file, particle_ID)
				particle_file_exist =  glob.glob(filein2)
				if not particle_file_exist:
					print 'File: ', filein, ' does not exist!'
					print 'Nor does file: ', filein2
					continue
				filein = filein2
		else :
			continue
		print "Plotting particle", j+1, 'of', len(particle_list)
#		pf = yt.load(infofile)
#		current_time = pf.current_time
		# If the output file in question exists, we then unpack
		rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, rhobin, mdotbin, norm, angXbin, angYbin, angZbin, vphi_magbin, sum_velocities, particleMass, creation_time, current_time, vrms_r_bin, vrms_l_bin, vrms_theta_bin, vrms_phi_bin = np.loadtxt(filein, unpack=True)

		creation_time = creation_time[0]
		current_time = current_time[0]
		if creation_time != 0.0:
			partAge = (current_time - creation_time) / yr / 1.e3
			print 'particle age in thousands of years', partAge
		else :
			partAge = 0.0
		# convert Particle Mass to cgs, mTbin is already in cgs
		particleMass = particleMass[0] * Msun
		print particleMass

		# Calculate the derivative of the density
		# And calculate the disk radius
		# Old method of finding disk radius
		#Density_derivative_profile(rbin, rhobin)
		#diskMass = mTbin[disk_radius_index] - particleMass
		#print 'Disk Mass is: ', diskMass, 'g or ', diskMass/Msun, 'M_sun'
		# New Method of finding disk radius.
		Disk_radius_by_vphi(rbin, vrbin, vrmsbin, vphi_magbin)
		diskMass = mTbin[disk_radius_index] - particleMass
		print 'Disk Mass is: ', diskMass, 'g or ', diskMass/Msun, 'M_sun'
		print 'Star particle Mass is: ', particleMass, 'g or ', particleMass/Msun, 'M_sun'
		try:
			print 'Disk/Star Mass Ratio is: ', diskMass / particleMass
		except ValueError:
			print 'No stellar mass so dividing by zero'
		disk_particle_mass = particleMass + diskMass
		r_star = Find_r_star(rbin, mTbin, particleMass)
		print 'R_* is: ', r_star, 'pcs or ', r_star * parsec, 'cms'
		MC15_velocity, py_mdot = Additional_Calculations()
		Thermal_P_ratio, Turb_P_ratio, Rot_P_ratio = Pressure_calculation()
		# If we're looking to get an average profile,
		if (args.avg):
			if len(plotting_list) > 1:
				print "please avg over one quantity at a time for now."
				sys.exit()
			particleMass_solar = particleMass/Msun
			lower_bound = center_mass_4_avg - deltaMass
			upper_bound = center_mass_4_avg + deltaMass
			if lower_bound <= particleMass/Msun <= upper_bound:
				print particleMass_solar
				savefig = False
				#for i in range(len(plotting_list)): # For now remove looping through all options.
				particle_count = particle_count + 1
				Setup_plot_window()
#				fileout = "avg_{0}_{1}Msun_{2}_{3}.{4}".format(plotting_list[i], center_mass_4_avg, compare_file, timestep, output_format)
#				if 'velocity' == plotting_list[i] :
				if( 'velocity' in plotting_list) :

					tot_vr, avg_vr = obtain_avg_value(tot_vr, particle_count, vrbin)
					print 'step', timestep
					print 'end', args.end-2
					if timestep >= (args.end-2):
						fileout = "avg_{0}_{1}Msun_{2}_{3}.{4}".format('velocity', center_mass_4_avg, compare_file, timestep, output_format)
						Velocity_plot(rbin, avg_vr, vrmsbin, vmagbin, vKbin, vphi_magbin, sum_velocities, MC15_velocity)
						savefig = True
					print len(rbin), len(vrbin), len(sigma)
					#print sigma[0]
#					plt.errorbar(rbin, vrbin, xerr, sigma[0])

#					print 'total', type(tot_vr), 'avg', type(avg_vr)
#					print tot_vr
#					print 'avg', avg_vr
#					
#				elif 'turbulent' == plotting_list[i] :
#					tot_vrms, avg_vrms = obtain_avg_value(tot_vrms, particle_count, vrmsbin)
#					Turbulent_velocity_plot(rbin, avg_vrms, vrms_r_bin, vrms_l_bin)
#					savefig = True
#				elif 'density' == plotting_list[i] :
#					tot_density, avg_density = obtain_avg_value(tot_density, particle_count, rhobin)
#					Density_plot(rbin, avg_density)
#					savefig = True
#				elif 'mdot' == plotting_list[i] :
#					#						py_mdot = -4.0*pi*rhobin*vrbin*(parsec*rbin)**2*yr/Msun) 
#					tot_mdot, avg_mdot = obtain_avg_value(tot_mdot, particle_count, py_mdot)
#					Mdot_plot(rbin, avg_mdot, mdotbin)
#					savefig = True
				if (savefig):
					png_save(fileout)
					savefig = False
			#If we are averaging we do not also want to 
			# write out all the plots, so go to next timestep.
			continue

		print 'Dispersing to plot.'
		# For the particle and Timestep we are currently in, loop and do all requested plots.
		for i in range(len(plotting_list)):
			print plotting_list[i]
			Setup_plot_window()
			fileout = "{0}_ID{1}_{2}_{3}.{4}".format(plotting_list[i], particle_ID, compare_file, timestep, output_format)
			if 'velocity' == plotting_list[i] :
				Velocity_plot(rbin, vrbin, vrmsbin, vmagbin, vKbin, vphi_magbin, sum_velocities, MC15_velocity)
			elif 'turbulent' == plotting_list[i] :
				Turbulent_velocity_plot(rbin, vrmsbin, vrms_r_bin, vrms_l_bin)
			elif 'toomreQ' == plotting_list[i] :
				Toomre_Q(rbin, rhobin, vrmsbin, vphi_magbin)
			elif 'veldens' == plotting_list[i] :
				Velocity_vs_Density_plot(rhobin, vrbin, vrmsbin, vKbin, vphi_magbin, particleMass, MC15_velocity)
			elif 'density' == plotting_list[i] :
				Density_plot(rbin, rhobin)
			elif 'angmomentum' == plotting_list[i] :
				Ang_moment_plot(rbin, vphi_magbin)
			elif 'pressure' == plotting_list[i] :
				Pressure_Gravity_ratio(rbin, rhobin, vrmsbin, vphi_magbin, Thermal_P_ratio, Turb_P_ratio, Rot_P_ratio)
			elif 'mass' == plotting_list[i] :
				Total_mass_plot(rbin, mTbin)
			elif 'mdot' == plotting_list[i] :
				Mdot_plot(rbin, py_mdot, mdotbin)
			png_save(fileout)
