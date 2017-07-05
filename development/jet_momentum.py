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
import yt
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


def standardize_File_number(input_File_number):
	output_File_number = "%05d"%input_File_number
	if len(output_File_number) >= 6:
		print 'This system assumes less than one million output files.'
		sys.exit()
	return output_File_number


############################################
############################################
############################################
parser = argparse.ArgumentParser(description = " Requires: start number, end number, stride_length, file_name, which_plots")
parser.add_argument('start', metavar='N1', type=int, help ='Start value for hdf5 files')
parser.add_argument('end', metavar='N2', type=int, help='End Value for hdf5 files, note runs until End-1')
parser.add_argument('step', metavar='N3', type=int, help='Stride length')
args = parser.parse_args()
global Init_matplotParams
Init_matplotParams = False
global Init_parts
Init_parts = False
file_prefix = 'rad_profile'
######### Constants #######
mp = 1.6e-24
pi = np.pi
parsec = 3e18
yr = 3.15 * 1e7
Msun = 2e33
G = 6.67e-8
c_s = 2.64e4
scale_t      =  0.387201000000000E+04

total_avg_v = []
time = []
for timestep in range(args.start,args.end,args.step) :   
	# For each timestep, check to see if the info and sink files exist.
	File_number = standardize_File_number(timestep)
	infofile = "output_{0}/info_{0}.txt".format(File_number)
	sinkfile = "output_{0}/sink_{0}.info".format(File_number)
	sinkcsvfile = "output_{0}/sink_{0}.csv".format(File_number)

	if( not os.path.isfile( sinkfile)):# or ( not os.path.isfile( infofile)):
		#print "Either the sink or info file for timestep", File_number, "is missing."
		print "The sink file for timestep", File_number, "is missing."
		continue #return
	try :
		particle_list, part_masses, part_pjet = np.loadtxt( sinkfile, usecols=[0,1,12], unpack=True, skiprows=3, comments="=")
		#pjet is loaded in units of Solar masses * km/s
		Particle_ID_list, partMass, xstar, ystar, zstar, \
		    vxstar, vystar, vzstar, lxstar, lystar, lzstar, \
		    Particle_Age, dMBhoverdt_star = np.loadtxt( sinkcsvfile, unpack=True, skiprows=0, delimiter=',', comments="=")
	except UserWarning :
		print 'Sink file is empty'
		print 'The sinkfile is empty and you have not requested to plot max density location.'
		#print "If you want to plot a particle before its formation, use the flag --noparticle."
		continue #return
	pf = yt.load(infofile) #Current time is loaded in COde
	current_time = np.float64(pf.current_time) / yr
	if not Init_parts:
		Init_parts = True
		print 'setting'
		start_time = current_time
	Particle_Age = Particle_Age * scale_t
	print File_number
	print current_time
	print start_time
	print current_time - start_time

	try:
		timestep_particles = len(particle_list)
	except TypeError:
		timestep_particles = 1
		part_pjet = np.array([part_pjet])
		part_masses = np.array([part_masses])
	time_step_m = 0
	time_step_p = 0
	print part_masses
	print part_pjet
        for index in range(timestep_particles):
            stellar_mass = part_masses[index]
            pjet = part_pjet[index]
            jet_mass = 0.5 * stellar_mass # solar masses.
#	    vjet = pjet / jet_mass # km/s

	    time_step_m = time_step_m + jet_mass
	    time_step_p = time_step_p + pjet
#	    time_step_v = time_step_v + vjet
	avg_v = time_step_p / time_step_m # mass weighted velocity.
#	print current_time / yr
#	print start_time / yr
#	print (current_time - start_time) /yr
	time.append(current_time - start_time)
	total_avg_v.append(avg_v)


np.savetxt('jet_momentum.txt', zip(time, total_avg_v), fmt="%15.9E")

plt.loglog(time, total_avg_v, 'b-', label='$u_r$')
#plt.ylim(1e-1, 1e2)
#plt.xlim(1e0, 1e6)
plt.ylabel(r'${\rm Jet \, velocity}$ $({\rm km\, s}^{-1})$', fontsize = 25)
plt.xlabel(r'$t-t_*$ $({\rm yr})$', fontsize=25)
plt.savefig('times_test.pdf')
