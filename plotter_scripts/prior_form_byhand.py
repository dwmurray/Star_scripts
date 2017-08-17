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

#####Constants##########
mp = 1.6e-24
pi = np.pi
parsec = 3e18
yr = pi * 1e7
Msun = 2e33
G = 6.67e-8
c_s = 2.64e4


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



def Update_avg_value_list(filein, particle_ID, File_number):
	""" This function determines if the particle is already accounted for, so as to avoid double counting. It does so by looking at the Particle ID that the sink particle is assigned in the simulation run.
	If the particle is already counted, but is closer to the value of mass that the user asked for, it will update the list with the closer value's information and remove the older value.
	"""
	global Particles_Used
	Particle_in_List = False
	#load in everything here, again, mostly to not pass in everything through the arguement list
	rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, rhobin, mdotbin, norm, angXbin, angYbin, angZbin, vphi_magbin, sum_velocities, particleMass, creation_time, current_time, vrms_r_bin, vrms_l_bin, vrms_theta_bin, vrms_phi_bin = np.loadtxt(filein, unpack=True)
	vrbin = np.nan_to_num(vrbin)
#	print vrbin
#	print ""
	vrbin = np.clip(vrbin,a_min=None,a_max=0)
#	print vrbin

	vrmsbin = np.nan_to_num(vrmsbin)
	vphi_magbin = np.nan_to_num(vphi_magbin)
	vKbin = np.nan_to_num(vKbin)
	rhobin = np.nan_to_num(rhobin)
	if Particles_Used == []:
		Particles_Used.append([particle_ID, File_number, deltaM, particleMass, rbin, vrbin, vrmsbin, vKbin, vmagbin, mTbin, rhobin, mdotbin, vphi_magbin, sum_velocities, vrms_r_bin, vrms_l_bin, vrms_theta_bin, vrms_phi_bin])
	else:
		for iter in range(len(Particles_Used)):
			sublist = Particles_Used[iter]#sublist is [particle_ID(int), File_number(str), deltaM(float64), particleMass, rbin, vrbin(np.ndarray), etc]
			if sublist[0] == particle_ID:
				Particle_in_List = True
				if deltaM < sublist[2]:
					Particles_Used[iter] = [particle_ID, File_number, deltaM, particleMass, rbin, vrbin, vrmsbin, vKbin, vmagbin, mTbin, rhobin, mdotbin, vphi_magbin, sum_velocities, vrms_r_bin, vrms_l_bin, vrms_theta_bin, vrms_phi_bin]
				break
		if not Particle_in_List:
			Particles_Used.append([particle_ID, File_number, deltaM, particleMass, rbin, vrbin, vrmsbin, vKbin, vmagbin, mTbin, rhobin, mdotbin, vphi_magbin, sum_velocities, vrms_r_bin, vrms_l_bin, vrms_theta_bin, vrms_phi_bin])
	print 'Current length of Part list : ', len(Particles_Used)



def Obtain_avg_value():
	""" This function calculates the average and standard deviation.
	"""
	global sigma
	avg_IDs = []
	avg_filenums = []
	particle_count = len(Particles_Used)
	for ob_ID in range(particle_count):
		if Particles_Used[ob_ID][0] not in avg_IDs:
			avg_IDs.append(Particles_Used[ob_ID][0])
			avg_filenums.append(Particles_Used[ob_ID][1])
	means = [avg_IDs, avg_filenums]
#	print 'particle count', particle_count
#	print 'size in table for each particle', len(Particles_Used[0])
	#Particles_Used.append([particle_ID[0], File_number[1], deltaM[2], particleMass[3], rbin[4], vrbin[5], vrmsbin[6], vKbin[7], vmagbin[8], mTbin[9], rhobin[10], mdotbin[11], vphi_magbin[12], sum_velocities[13], vrms_r_bin[14], vrms_l_bin[15], vrms_theta_bin[16], vrms_phi_bin[17]])
	for index_2_average in range(2,len(Particles_Used[0])): #start at 2 to avoid avging over the particle IDs and the timestamp file numbers.
#		print index_2_average
		tot_bin = 0. # reset total bin
		for iter in range(len(Particles_Used)):
#			print Particles_Used[iter][index_2_average]
			tot_bin = tot_bin + Particles_Used[iter][index_2_average]
#		print tot_bin
		mean = tot_bin/particle_count
#		print 'mean', mean
		means.append(mean)
#	sigma = np.sqrt((Values_4_std - mean)**2/particle_count)
#	print len(Values_4_std), len(sigma), particle_count
#	print 'Values_4_std', Values_4_std
#	print 'sigma', sigma
	return means

############################################
parser = argparse.ArgumentParser(description = " Requires: start number, end number, stride_length, file_name, which_plots")
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
parser.add_argument('--magnetic', action='store_true')
parser.add_argument('--Q', action='store_true')
parser.add_argument('--pjet', action='store_true')

parser.add_argument('--eighty', action='store_true')
parser.add_argument('--forty', action='store_true')

args = parser.parse_args()

plotting_list = []
accepted_plotting_list = ['velocity', 'turbulent', 'toomreQ', 'veldens', 'density', 'angmomentum', 'pressure', 'mass', 'mdot']#, 'magnetic']
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
	if ( args.magnetic):
		plotting_list.append('magnetic')
if plotting_list == []:
	print "You have not selected the plots to create."
	print "Possible plots include: ", accepted_plotting_list
	sys.exit()

global Init_matplotParams
Init_matplotParams = False

Particles_Used = []
#-40kyr
if args.forty :
	files = ['rad_profile_00120_shellsphere_1.out', 'rad_profile_00120_shellsphere_2.out', 'rad_profile_00120_shellsphere_3.out', \
			 'rad_profile_00122_shellsphere_4.out', 'rad_profile_00123_shellsphere_5.out', 'rad_profile_00123_shellsphere_6.out', \
			 'rad_profile_00123_shellsphere_7.out', 'rad_profile_00124_shellsphere_8.out', 'rad_profile_00124_shellsphere_9.out', \
			 'rad_profile_00124_shellsphere_10.out', 'rad_profile_00124_shellsphere_11.out', 'rad_profile_00125_shellsphere_12.out', \
			 'rad_profile_00128_shellsphere_13.out', 'rad_profile_00128_shellsphere_14.out', 'rad_profile_00129_shellsphere_15.out']


#-80kyr
if args.eighty :
	files = ['rad_profile_00119_shellsphere_1.out', 'rad_profile_00119_shellsphere_2.out', 'rad_profile_00119_shellsphere_3.out', \
			 'rad_profile_00121_shellsphere_4.out', 'rad_profile_00122_shellsphere_5.out', 'rad_profile_00122_shellsphere_6.out', \
			 'rad_profile_00122_shellsphere_7.out', 'rad_profile_00123_shellsphere_8.out', 'rad_profile_00123_shellsphere_9.out', \
			 'rad_profile_00123_shellsphere_10.out', 'rad_profile_00123_shellsphere_11.out', 'rad_profile_00124_shellsphere_12.out', \
			 'rad_profile_00127_shellsphere_13.out', 'rad_profile_00127_shellsphere_14.out', 'rad_profile_00128_shellsphere_15.out']

#files = ['rad_profile_00119_shellsphere_1.out', 'rad_profile_00119_shellsphere_2.out', 'rad_profile_00119_shellsphere_3.out', 'rad_profile_00121_shellsphere_4.out', 'rad_profile_00122_shellsphere_5.out', 'rad_profile_00122_shellsphere_6.out']
for file in files:
	particle_ID = 42
	File_number = 42
	deltaM = 42
	rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, rhobin, \
	    mdotbin, norm, angXbin, angYbin, angZbin, vphi_magbin, sum_velocities, \
	    particleMass, creation_time, current_time, vrms_r_bin, vrms_l_bin, vrms_theta_bin, vrms_phi_bin = np.loadtxt(file, unpack=True)
	vrbin = np.nan_to_num(vrbin)
	vrbin = np.clip(vrbin,a_min=None,a_max=0)
	vrmsbin = np.nan_to_num(vrmsbin)
	vphi_magbin = np.nan_to_num(vphi_magbin)
	vKbin = np.nan_to_num(vKbin)
	vmagbin = np.nan_to_num(vmagbin)
	rhobin = np.nan_to_num(rhobin)
	mdotbin = np.nan_to_num(mdotbin)
	sum_velocities = np.nan_to_num(sum_velocities)
	vrms_l_bin = np.nan_to_num(vrms_l_bin)
	vrms_r_bin = np.nan_to_num(vrms_r_bin)
	vrms_theta_bin = np.nan_to_num(vrms_theta_bin)
	vrms_phi_bin = np.nan_to_num(vrms_phi_bin)
	Particles_Used.append([particle_ID, File_number, deltaM, particleMass, rbin, vrbin, vrmsbin, vKbin, vmagbin, mTbin, rhobin, mdotbin, vphi_magbin, sum_velocities, vrms_r_bin, vrms_l_bin, vrms_theta_bin, vrms_phi_bin])

means = Obtain_avg_value()
#	print means
print -4.0*pi*means[10]*means[5]*(parsec*means[4])**2*yr/Msun
compare_file = 'shellsphere'
output_format = 'png'	
for i in range(len(plotting_list)):
	if args.eighty :
		savefileout = "avg_{0}_by_hand_{1}_80kyr.txt".format(plotting_list[i], compare_file, output_format)
	elif args.forty :
		savefileout = "avg_{0}_by_hand_{1}_40kyr.txt".format(plotting_list[i], compare_file, output_format)
	if 'velocity' == plotting_list[i] :
		np.savetxt(savefileout, zip(means[4], means[5], means[6], means[8], means[7], means[12]), fmt="%15.9E")
	elif 'density' == plotting_list[i] :
		np.savetxt(savefileout, zip(means[4], means[10]), fmt="%15.9E")
	elif 'mdot' == plotting_list[i] :
		np.savetxt(savefileout, zip(means[4], -4.0*pi*means[10]*means[5]*(parsec*means[4])**2*yr/Msun), fmt="%15.9E")


