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
#rcParams['ytick.labelsize']=25

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

def density_original_plot(rbin1, mTbin1, rhobin1, partAge1, rbin2, mTbin2, rhobin2, partAge2):
	if partAge1 < 0.0:
		particle_1_label = str(int(abs(partAge1))) + ' yr prior to formation'
	else:
		particle_1_label = str(int(partAge1)) + ' yr after formation'

	if partAge2 < 0.0:
		particle_2_label = str(int(abs(partAge2))) + ' yr prior to formation'
	else:
		particle_2_label = str(int(partAge2)) + ' yr after to formation'

	pl.clf()	
	pl.loglog(rbin1, rhobin1/mp, 'b.--', label = particle_1_label)
	pl.loglog(rbin2, rhobin2/mp, 'g-', label = particle_2_label)
	pl.legend(loc=0, fontsize='22', frameon=False)
	pl.xlim(3e-3, 3e0)
	pl.ylim(1e0, 1e10)
	#pl.xticks(fontsize= 16)
	#pl.yticks(fontsize= 16)
	pl.tick_params('both',length=8, width=1, which='minor')
	pl.tick_params('both',length=10, width=1.5, which='major')
	pl.gcf().subplots_adjust(bottom=0.15)
	pl.ylabel('$N$ ($cm^{-3}$)', fontsize = 16)
	pl.xlabel('$r$ ($pc$)', fontsize = 16)
	
def density_test_plot(rbin1, mTbin1, rhobin1, partAge1, rbin2, mTbin2, rhobin2, partAge2):#(rbin, mTbin, rhobin):
	pl.clf()
	pl.rc('text', usetex=True)
	pl.rc('font', family='serif')

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
	pl.xticks(fontsize=30)
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
	host.loglog(rbin1, rhobin1/mp, label='Number density') 	
	host.loglog(rbin2, rhobin2/mp, label='Number density') 	
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
	#pl.title('Particle {0} {1:03d} kyr after formation'.format(particle_name, int(partAge)))
	#pl.title('Run of Density {0:03d} kyr after formation'.format(int(partAge/1e3)))
	pl.rc('text', usetex=False)

def density_plot(rbin1, mTbin1, rhobin1, partAge1, rbin2, mTbin2, rhobin2, partAge2):
	if partAge1 < 0.0:
		particle_1_label = str(int(abs(partAge1))) + ' yr prior to formation'
	else:
		particle_1_label = str(int(partAge1)) + ' yr after formation'

	if partAge2 < 0.0:
		particle_2_label = str(int(abs(partAge2))) + ' yr prior to formation'
	else:
		particle_2_label = str(int(partAge2)) + ' yr after formation'

	pl.clf()
	pl.rc('text', usetex=True)
	pl.rc('font', family='serif')
	host = host_subplot(111, axes_class=AA.Axes)
	par1 = host.twinx()
	Ndensity_Y_min = 1e1
	Ndensity_Y_max = 1e8
	host.set_xlim(3e-3, 5e0)
	host.set_ylim(Ndensity_Y_min, Ndensity_Y_max)
	Mdensity_Y_min = Ndensity_Y_min * 2. * mp
	Mdensity_Y_max = Ndensity_Y_max * 2. * mp
	par1.set_ylim(Mdensity_Y_min, Mdensity_Y_max)
	par1.set_yscale('log')
	pl.gcf().subplots_adjust(bottom=0.15)
	host.set_ylabel('$n$ $({\\rm cm}^{-3})$', fontsize = 28)
	host.set_xlabel('$r$ $({\\rm pc})$', fontsize = 28)
	par1.set_ylabel('$\\rho$ $({\\rm g\\, cm}^{-3})$', fontsize = 28)
	host.axis["left"].label.set_fontsize(25)
	host.axis["bottom"].label.set_fontsize(25)
	par1.axis["right"].label.set_fontsize(25)
	host.loglog(rbin1, rhobin1/mp, 'b.--', label = particle_1_label)
	host.loglog(rbin2, rhobin2/mp, 'g-', label = particle_2_label)
	pl.legend(loc=0, fontsize='20', frameon=False)
	pl.rc('text', usetex=False)

def png_save(png_fileout):
	print 'saving file: ', png_fileout
	#pl.draw()
	pl.savefig(png_fileout)

parser = argparse.ArgumentParser(description = " Requires: start number, end number, stride_length, file_name, which_plots")

parser.add_argument('first', metavar='N1', type=int)
parser.add_argument('second', metavar='N2', type=int)
parser.add_argument('ParticleID', metavar='N3', type=int)

args = parser.parse_args()
file1 = args.first
file2 = args.second
particle_number = args.ParticleID
#particle_number = particle_number + '.0'

mp = 1.6e-24
pi = np.pi
parsec = 3e18
sec_in_yr = np.pi* 1e7
Msun = 2e33
G = 6.67e-8
c_s = 2.64e4
Gravity_Turned_On = 6.317465983e14

compare_files = 'shellsphere'
file_prefix = 'rad_profile'
quad = os.getcwd()[-5:]

print particle_number
print type(particle_number)
# The file looks like this:
#'{file_prefix}{framestep}_{compare_files}_{particle_number}.out'
# i.e. rad_profile_0218_part_000.out
#file1_in = '{0}_{1}_{2:04d}_{3}_*.out'.format(file_prefix, quad, file1, compare_files)
#Flash:
file_load_1 = '{0}_{1}_{2:04d}_{3}_{4}.out'.format(file_prefix, quad, file1, compare_files, particle_number)
file_load_2 = '{0}_{1}_{2:04d}_{3}_{4}.out'.format(file_prefix, quad, file2, compare_files, particle_number)
# RAMSES:
#file_load_1 = '{0}_{1:04d}_{2}_{3}.out'.format(file_prefix, file1, compare_files, particle_number)
file_load_1 = '{0}_{1:04d}_{2}_{3}.0.out'.format(file_prefix, file1, compare_files, particle_number)
#file_load_2 = '{0}_{1:04d}_{2}_{3}.out'.format(file_prefix, file2, compare_files, particle_number)
file_load_2 = '{0}_{1:04d}_{2}_{3}.0.out'.format(file_prefix, file2, compare_files, particle_number)

file_exist1 = glob.glob(file_load_1)
if file_exist1:
	rbin1, vrbin1, vrmsbin1, vrmsnbin1, vKbin1, vmagbin1, vmagnbin1, mTbin1, rhobin1, mdotbin1, norm1, angXbin1, angYbin1, angZbin1, vphi_magbin1, sum_of_velocities1, partMass1, part_creation_time1, current_time1, vrms_r_bin1, vrms_l_bin1, vrms_theta_bin1, vrms_phi_bin1  = np.loadtxt(file_load_1, unpack=True)
	partCreationtime1 = part_creation_time1[0]
	print partCreationtime1
#	print rbin1, rhobin1
	current_time1 = current_time1[0]
	partAge1 = current_time1 - partCreationtime1
	partAge1 = partAge1 / sec_in_yr
else:
	print 'File', file_load_1, " does not exist, exiting gracefully"
	print 'File doesn"t exist, exiting gracefully'
	sys.exit()

file_exist2 = glob.glob(file_load_2)
if file_exist2:
	rbin2, vrbin2, vrmsbin2, vrmsnbin2, vKbin2, vmagbin2, vmagnbin2, mTbin2, rhobin2, mdotbin2, norm2, angXbin2, angYbin2, angZbin2, vphi_magbin2, sum_of_velocities2, partMass2, part_creation_time2, current_time2, vrms_r_bin2, vrms_l_bin2, vrms_theta_bin2, vrms_phi_bin2  = np.loadtxt(file_load_2, unpack=True)
	partCreationtime2 = part_creation_time2[0]
	print partCreationtime2
	current_time2 = current_time2[0]
	partAge2 = current_time2 - partCreationtime2
	partAge2 = partAge2 / sec_in_yr
else:
	print 'File', file_load_2, ' doesn"t exist, exiting gracefully'
	sys.exit()

if current_time2 >= current_time1:
	print 'time since gravity (kyr): ', (current_time2 - Gravity_Turned_On) / (3.14e10)
elif current_time1 >= current_time2:
	print 'time since gravity (kyr): ', (current_time1 - Gravity_Turned_On) / (3.14e10)

png_fileout = "double_density_{0}_{1:04d}_{2:04d}_{3}.png".format(quad, file1, file2, args.ParticleID)
density_plot(rbin1, mTbin1, rhobin1, partAge1, rbin2, mTbin2, rhobin2, partAge2)
png_save(png_fileout)
