import matplotlib
matplotlib.use("Agg")
from matplotlib import rc
#rc("text", usetex=True)
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
#

import h5py
import numpy
import os
import sys
import glob
import matplotlib.pyplot as pl

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

def run_quad( quad) :
	times = []
	masses = []
	number_of_stars = []
	part_prefix = 'BB_hdf5_part_'

	for fileIndex in range( 1, 1000) : 
		filename = "{0}/{1}{2:04d}".format(quad,part_prefix, fileIndex) 
		#filename = "q2_{0}/{1}{2:04d}".format(quad,part_prefix, fileIndex) 
		#filename = "jeans_4096_scinet/{0}{1:04d}".format(part_prefix, fileIndex)
		#filename = "jeans_padoan_feder/{0}{1:04d}".format(part_prefix, fileIndex)
		#filename = "jeans_1024/{0}{1:04d}".format(part_prefix, fileIndex)
		file_exist = glob.glob(filename)
		if not file_exist:
			#print 'File: '+ filename + '  does not exist!'
			continue
	
		particleFile = File_opener(filename)

		if( particleFile != None) : 
			current_time = particleFile['real scalars'][0][1]
			
			partIndices = getValues(particleFile, "tag").astype("int")
			#print partIndices
			xps = getValues(particleFile, "posx")
			yps = getValues(particleFile, "posy")
			zps = getValues(particleFile, "posz")
			vx = getValues(particleFile, "velx")
			vy = getValues(particleFile, "vely")
			vz = getValues(particleFile, "velz")
			
			mass_part_solar = getValues(particleFile, "mass")/Msun
			partCreateTime = getValues(particleFile, "creation_time")
		
			times.append(current_time)
			masses.append(mass_part_solar.sum())
			particleFile.close()
			number_of_stars.append(len(mass_part_solar))
			#number_of_stars =len(mass_part_solar)
			sortList = numpy.argsort( mass_part_solar)[::-1]
			print " "
			#print "Masses:" + str(zip(partIndices[sortList],mass_part_solar[sortList]))
			print "Masses:" + str(mass_part_solar[sortList])
#			print number_of_stars
	return numpy.array(times), numpy.array(masses), number_of_stars

# Constants #
mp = 1.6e-24
pi = numpy.pi
parsec = 3e18
sec_in_yr = numpy.pi* 1e7
Msun = 2e33
G = 6.67e-8
c_s = 2.64e4

qtimes = []
qmasses = []
totnumber_stars = []
#for quad in range(1, 6) : 
#	qtime, qmass = run_quad( quad)
#	qtimes.append(qtime)
#	qmasses.append(qmass)
#	print quad
#	print qtime[-1]
#
#quad_list = [1,2,3,4,5]
#quad_list = [2,3]
quad_list = [2]
for j in range(0, 1) : 
	quad = quad_list[j]
	print quad
	qtime, qmass, number_stars = run_quad( quad)
	qtimes.append(qtime)
	qmasses.append(qmass)
	totnumber_stars.append(number_stars)
	print 'time of 1st particle', qtime[0]
	print 'last time in octant', qtime[-1]
# find the start, end time
start_time = 1e30
end_time = 1e30
for qtime  in qtimes : 
	if( qtime[0] < start_time) :
		start_time = qtime[0]
		print 'particle start time of the box', start_time
	if( qtime[-1] < end_time) : 
		end_time = qtime[-1]
		print 'The earliest end time', end_time



bins = 1001
times = numpy.arange(0, bins, 1) * (end_time-start_time)/bins + start_time
#times = numpy.arange( start_time, end_time, (end_time - start_time)/bins)
print times.size, start_time, end_time, end_time - start_time
import scipy.interpolate as si
totalMass = numpy.zeros( bins) 
totalnumber_stars = numpy.zeros(bins)
for qtime, qmass, number_stars  in zip( qtimes, qmasses, totnumber_stars) :
	
	func = si.interp1d(qtime, qmass, bounds_error = False, fill_value = 0e0, kind=3) 
	totalMass = totalMass +  func(times)
	func2 = si.interp1d(qtime, number_stars, bounds_error = False, fill_value = 0e0, kind=3) 
	totalnumber_stars = totalnumber_stars  + func2(times)

for i in qmass:
	print i
	#print totalMass[i]

times -= start_time
times_yr = times / (numpy.pi*1.e7)
times_hkyr = times_yr /(1.e5)


#print len(times_yr), len(totalnumber_stars)
#print totalnumber_stars
#pl.loglog( times_yr, totalMass) 
pl.loglog( (qtime -start_time)/3.15e7, qmass) 
#pl.loglog( times_yr, totalnumber_stars) 
#pl.loglog( times_hkyr, (times_hkyr/1e12)**2) 
pl.loglog( numpy.arange(2, 5, 0.1)*1e4, 3.*(numpy.arange(2,5,0.1))**2, ls="--") 
#pl.loglog( numpy.arange(2, 5, 0.1)*1e4, 4.*(numpy.arange(2,5,0.1))**1, ls="-.-") 
#pl.ylim(1.2e0, 3e2)
#pl.xlim(1.e4, 6e5)
#pl.ylim(1.2e-1, 3e2)
#pl.xlim(1.e3, 6e5)
#pl.xticks(fontsize=25)
#pl.yticks(fontsize=25)
#pl.tick_params('both',length=8, width=1, which='minor')
#pl.tick_params('both',length=10, width=1.5, which='major')
#pl.gcf().subplots_adjust(bottom=0.15)
pl.ylabel('$M_*$ $(M_\odot)$', fontsize = 25)
pl.xlabel( "$t-t_*$ $({\\rm yrs})$", fontsize = 25)
pl.savefig("Total_mass_over_time.pdf")

for t,m in zip((qtime-start_time)/3.15e7, qmass) :
	print "{0} {1}".format(t,m)
#pl.savefig("Total_mass_over_time.png")


#pl.clf()
#pl.loglog( times_yr, totnumber_stars)
#pl.loglog( numpy.arange(2, 5, 0.1)*1e4, 3.*(numpy.arange(2,5,0.1))**2, ls="--") 
#pl.loglog( numpy.arange(2, 5, 0.1)*1e4, 4.*(numpy.arange(2,5,0.1))**1, ls="-.-") 
#pl.ylabel('$N_*$ $$', fontsize = 25)
#pl.xlabel( "$t-t_*$ $({\\rm yrs})$", fontsize = 25)
#pl.savefig("Total_number_over_time.pdf")
