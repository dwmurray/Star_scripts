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


def obtain_avg_value(avg_bin, part_count_avg_over, quantity_to_avg_bin):
	part_count_avg_over = part_count_avg_over + 1
	avg_bin = avg_bin + quantity_to_avg_bin
	return avg_bin, part_count_avg_over, avg_bin/part_count_avg_over

	
def density_plot(rbin, mTbin, rhobin, avg_bin, particle_count):
	pass
#return avg_bin, particle_count

def png_save(fileout):
	print 'saving file: ', fileout
	pl.savefig(fileout)

#
#def getValues( particleFile, value) :
#        names = np.char.strip(particleFile["particle names"].value)
#        starParticles = particleFile["tracer particles"].value
#        ind = (np.char.strip(names) == value)[:,0]
#        return starParticles[:, ind].flatten()

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
	return partIndices

def Obtain_particles(sinkfile) :
        if( not os.path.isfile( sinkfile)) :
                return
        try :
                ID = np.loadtxt( sinkfile, usecols=[0], unpack=True, skiprows=3, comments="=")
		print ID
        except ValueError :
                print 'Empty sinkfile'
	if withNoParticle:
		ID = str(withParticleIDValue)
		print 'We call to obtain particles'
        return ID

parser = argparse.ArgumentParser(description = " Requires: start number, end number, stride_length, file_name, which_plots")

parser.add_argument('start', metavar='N1', type=int)
parser.add_argument('end', metavar='N2', type=int)
parser.add_argument('step', metavar='N3', type=int)
parser.add_argument('ParticleID', metavar='N4', type=int, nargs='?', default=42, help='Particle ID you want to reduce.')
parser.add_argument('--shellsphere', action='store_true')
parser.add_argument('--allparticles', action='store_true')
parser.add_argument('--allplots', action='store_true')
parser.add_argument('--disk', action='store_true')
parser.add_argument('--Q', action='store_true')
parser.add_argument('--pressure', action='store_true')
parser.add_argument('--pdf', action='store_true')

args = parser.parse_args()

withAllParticles = args.allparticles
withShellSphere = args.shellsphere
withParticleIDValue = args.ParticleID
withPDF = args.pdf
withNoParticle = False
mp = 1.6e-24
pi = np.pi
parsec = 3e18
sec_in_yr = np.pi* 1e7
Msun = 2e33
G = 6.67e-8
c_s = 2.64e4

Gravity_Turned_On = 6.317465983e14

# These have to do with calculating the sphere of influence
chosen_ratio_number = 3.0
particleMass = 0.0
particle_count = 0
avg_bin = np.zeros(70)

compare_files = 'shellsphere_'

if (withPDF):
	output_format = 'pdf'
else:
	output_format = 'png'

file_prefix = 'rad_profile'
quad = os.getcwd()[-2:]
#quad = 'hd32'
quad = 'hd8'
#quad = 'hd4'

for framestep in range(args.start,args.end,args.step) :   
	sinkfile = "sink_{0:05d}.info".format(framestep)
	print sinkfile
	particle_list = Obtain_particles(sinkfile)
	if not withNoParticle:
		if particle_list.size == 1:
			print 'Particles have not been created yet'
			print 'Setting the particle list to include the specified particle id.'
			particle_list = np.array([withParticleIDValue])
			print particle_list
			print type(particle_list)
		else:
			print 'The total number of star particles in this timestep is: ', len(particle_list)
	for j in range(len(particle_list)):
		particle_number = int(particle_list[j])
		if particle_number == withParticleIDValue or (withAllParticles):
			print "Plotting particle", j+1, 'of', len(particle_list)
			filein = '{0}_{1:04d}_{2}{3}.0.out'.format(file_prefix, framestep, compare_files, particle_number)
			particle_file_exist =  glob.glob(filein)
			print filein
			if not particle_file_exist:
				filein = '{0}_{1:04d}_{2}{3}.out'.format(file_prefix, framestep, compare_files, particle_number)
				particle_file_exist =  glob.glob(filein)
				print filein
				if not particle_file_exist:
					print 'File: ', filein, ' does not exist!'
					continue
			# If the file in question exists, we then unpack
			rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, rhobin, mdotbin, norm, angXbin, angYbin, angZbin, vphi_magbin, sum_of_velocities, partMass, part_creation_time, current_time, vrms_r_bin, vrms_l_bin, vrms_theta_bin, vrms_phi_bin = np.loadtxt(filein, unpack=True)

			partMass_solar = partMass #save the particle mass in solar masses for the avg.
			partMass = partMass * Msun
			
			# If only one particle set =, if a list, pick which particle its on
			if isinstance(partMass, float) ==True:
				particleMass = partMass
			else:
				particleMass = partMass[j]
			#r_star = 0.0
			fileout = "density_{0}_{1:04d}_{2}avg.{4}".format(quad, framestep, compare_files, particle_number, output_format)
			print "particle_count sent to avg", particle_count

			print 'particle_count at this location:', particle_count
			print partMass_solar[0]
			wanted_mass = 1.0
			print wanted_mass - wanted_mass*0.5
			
			if 0.5 <= partMass_solar[0] <= 3.0: #wanted_mass + wanted_mass * 0.5:
				#print rhobin
				#print type(rhobin)
				#rhobin = np.array(rhobin)
				#np.nan_to_num(x)
				rhobin = np.nan_to_num(rhobin)
				#print rhobin
				avg_bin, particle_count, avg_rhobin = obtain_avg_value(avg_bin, particle_count, rhobin)

				print avg_rhobin/mp
				print particle_count
				print ""
				avg_dens_fileout = "avg_density_{0}_{1:04d}.out".format(quad, framestep)
				np.savetxt(avg_dens_fileout, zip(rbin, avg_rhobin/mp), fmt="%15.9E")
				print 'particle_count at this location:', particle_count
		else:
			print "The particle ID's don't match, proceeding to next particle."
			continue
