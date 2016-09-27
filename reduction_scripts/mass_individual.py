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


import h5py
import numpy
import os
import sys
import argparse
import glob
import matplotlib.pyplot as pl
import scipy.interpolate as si

def getValues( particleFile, value) :
        names = numpy.char.strip(particleFile["particle names"].value)
        starParticles = particleFile["tracer particles"].value
        ind = (numpy.char.strip(names) == value)[:,0]
        return starParticles[:, ind].flatten()

def File_opener(filename):
    #particleFile = h5py.File("{0}{1:04d}".format(part_prefix,fileIndex))
    particleFile = h5py.File(filename)
#    print particleFile.keys()
# get all the star particles first
    if("tracer particles" in particleFile.keys()) :
	    return particleFile
    else :
	    particleFile.close()
	    return None

def run_quad( quad, withParticleIDValue) :
	times = []
	masses = []
	part_prefix = 'BB_hdf5_part_'
	for fileIndex in range( 1, 1000) :
		filename = "quad{0}/{1}{2:04d}".format(quad,part_prefix, fileIndex) 
		file_exist = glob.glob(filename)
		if not file_exist:
			continue
		particleFile = File_opener(filename)
		if( particleFile != None) : 
			current_time = particleFile['real scalars'][0][1]
			partIndices = getValues(particleFile, "tag").astype("int")
			if withParticleIDValue in partIndices:
				pass
			else:
				#print 'Particle {0} not in this file'.format(withParticleIDValue)
				particleFile.close()
				continue
			xps = getValues(particleFile, "posx")
			yps = getValues(particleFile, "posy")
			zps = getValues(particleFile, "posz")
			vx = getValues(particleFile, "velx")
			vy = getValues(particleFile, "vely")
			vz = getValues(particleFile, "velz")
			mass_part_solar = getValues(particleFile, "mass")/Msun
			partCreateTime = getValues(particleFile, "creation_time")
			for j in range(len(partIndices)):
				if partIndices[j] == withParticleIDValue:
					times.append(current_time)
					masses.append(mass_part_solar[j])
			particleFile.close()

	return numpy.array(times), numpy.array(masses)

def get_particle_list( quad) :
	global last_particle_list_quad1, last_particle_list_quad2, last_particle_list_quad3
	global last_particle_list_quad4, last_particle_list_quad5
	part_prefix = 'BB_hdf5_part_'
	for fileIndex in range( 1, 1000) : 
		filename = "quad{0}/{1}{2:04d}".format(quad,part_prefix, fileIndex) 
		file_exist = glob.glob(filename)
		if not file_exist:
			continue
	
		particleFile = File_opener(filename)
		if( particleFile != None) : 
			partIndices = getValues(particleFile, "tag").astype("int")
			if quad == 1:
				last_particle_list_quad1 = partIndices
			elif quad == 2:
				last_particle_list_quad2 = partIndices
			elif quad == 3:
				last_particle_list_quad3 = partIndices
			elif quad == 4:
				last_particle_list_quad4 = partIndices
			elif quad == 5:
				last_particle_list_quad5 = partIndices
			particleFile.close()

def plotting(qtime, qmass, start_time, end_time, withParticleIDValue):
	pl.clf()
	#bins = 1001
	#times = numpy.arange( start_time, end_time, (end_time - start_time)/bins)
	#print times.size, start_time, end_time, end_time - start_time

	#totalMass = numpy.zeros( bins) 
	#for qtime, qmass in zip( qtimes, qmasses) :
	#func = si.interp1d(qtime, qmass, bounds_error = False, fill_value = 0e0, kind=3) 
	#totalMass = totalMass +  func(times)
	#print totalMass
	times = qtime
	times -= start_time
	times_yr = times / (numpy.pi*1.e7)
	#times_hkyr = times_yr /(1.e5)
	#pl.loglog( times_yr, totalMass) 
	pl.loglog(times_yr, qmass) 
	#pl.loglog( times_hkyr, (times_hkyr/1e12)**2) 
	pl.loglog( numpy.arange(2, 50, 0.1)*1e4, 4.*(numpy.arange(2,50,0.1))**2, ls="--") 
	pl.loglog( numpy.arange(2, 50, 0.1)*1e4, 4.*(numpy.arange(2,50,0.1)), ls="-.") 
	#pl.ylim(1.2e0, 3e2)
	#pl.xlim(1.e4, 6e5)
	pl.xticks(fontsize=25)
	pl.yticks(fontsize=25)
	#pl.tick_params('both',length=8, width=1, which='minor')
	#pl.tick_params('both',length=10, width=1.5, which='major')
	pl.gcf().subplots_adjust(bottom=0.15)
	pl.ylabel('$M_*$ ($M_\odot$)', fontsize = 25)
	pl.xlabel( "$t-t_*$ (${\\rm yrs}$)", fontsize = 25)
	#pl.savefig("Total_mass_over_time_all.pdf")
	pl.savefig("Total_mass_over_time_ID{0}.pdf".format(withParticleIDValue))



#parser = argparse.ArgumentParser(description = " Requires: start number, end number, stride_length, file_name, which_plots")
#parser.add_argument('ParticleID', metavar='N1', type=int, nargs='?', default=42, help='Particle ID you want to get m(t).')
#args = parser.parse_args()
#withParticleIDValue = args.ParticleID

# Constants #
mp = 1.6e-24
pi = numpy.pi
parsec = 3e18
sec_in_yr = numpy.pi* 1e7
Msun = 2e33
G = 6.67e-8
c_s = 2.64e4

last_particle_list_quad1 = []
last_particle_list_quad2 = []
last_particle_list_quad3 = []
last_particle_list_quad4 = []
last_particle_list_quad5 = []
for quad in range(1, 6) :
	get_particle_list(quad)

total_particle_list = []
total_particle_list.append(last_particle_list_quad1)
total_particle_list.append(last_particle_list_quad2)
total_particle_list.append(last_particle_list_quad3)
total_particle_list.append(last_particle_list_quad4)
total_particle_list.append(last_particle_list_quad5)
#sys.exit()
print len(last_particle_list_quad1)
print len(last_particle_list_quad2)
print len(last_particle_list_quad3)
print len(last_particle_list_quad4)
print len(last_particle_list_quad5)
print len(last_particle_list_quad1) + len(last_particle_list_quad2) + len(last_particle_list_quad3) + len(last_particle_list_quad4) +len(last_particle_list_quad5)
particle_counter = 0
for quad in range(1, 6):
	current_particle_list = total_particle_list[quad-1]
	print current_particle_list
	for current_particle in range(len(current_particle_list)):
		withParticleIDValue = current_particle_list[current_particle]
		print quad, withParticleIDValue
		qtime, qmass = run_quad(quad, withParticleIDValue)
		# find the start, end time
		start_time = qtime[0]
		end_time = qtime[-1]
		plotting(qtime, qmass, start_time, end_time, withParticleIDValue)
		particle_counter += 1

print particle_counter
