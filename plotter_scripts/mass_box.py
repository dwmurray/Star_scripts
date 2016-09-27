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
import yt
import sys
import glob
import matplotlib.pyplot as pl

def getValues( particleFile, value) :
        names = numpy.char.strip(particleFile["particle names"].value)
        starParticles = particleFile["tracer particles"].value
        ind = (numpy.char.strip(names) == value)[:,0]
        return starParticles[:, ind].flatten()

def Obtain_particles(file, sinkfile, fileIndex) :
	"""
	This function is used to extract the relevant information from the sink_###.info ASCI file produced by RAMSES.

	"""
        if( not os.path.isfile( sinkfile)) :
		#print os.path.isfile(sinkfile)
                return
        if( not os.path.isfile(file)) :
		#print os.path.isfile(file)
		return
        try :
                ID, mass = numpy.loadtxt( sinkfile, usecols=[0, 1], unpack=True, skiprows=3, comments="=")
		#print ID, mass
        except ValueError :
                print 'Empty sinkfile'
		ID = []
		mass = []
	pf = yt.load(file)
	current_time = float(pf.current_time)
	#print("Time unit: ", pf.time_unit)
	#current_time = float(pf.current_time.in_units('s'))
	#current_time = 0.0
        return ID, mass, current_time

def run_quad(dir_label):
	times = []
	masses = []
	number_of_stars = []
	prefix = 'output_'
	start_time = 0
	calc_time = 0
	print 'looking at all sink files in dir.', dir_label
	for fileIndex in range( 1, 500) : 
		file = prefix+"{0:05d}/info_{0:05d}.txt".format(fileIndex)
		sinkfile = prefix+"{0:05d}/sink_{0:05d}.info".format(fileIndex)
		if not glob.glob(file):
			continue
		if not glob.glob(sinkfile):
			continue
		print file, sinkfile
		ID, mass_part_solar, current_time = Obtain_particles(file, sinkfile, fileIndex)
		if ID == []:
			pass
		else :
			if start_time != 0:
				times.append(current_time)
				masses.append(mass_part_solar.sum())
				end_time = current_time
			if start_time == 0 :
				start_time = current_time
				print start_time, 'and mass is', mass_part_solar
				times.append(current_time)
				masses.append(mass_part_solar.sum())
			try:
				number_of_stars.append(len(mass_part_solar))
			except TypeError:
				number_of_stars = [1]
				#number_of_stars = len(mass_part_solar)
#				sortList = numpy.argsort( mass_part_solar)[::-1]
#				print " "
#				#print "Masses:" + str(zip(partIndices[sortList],mass_part_solar[sortList]))
#				print "Masses:" + str(mass_part_solar[sortList])
#				#print number_of_stars
#			else:
				#times.append(current_time)
				#number_of_stars.append(len(mass_part_solar))
				#number_of_stars =len(mass_part_solar)
				#sortList = numpy.argsort( mass_part_solar)[::-1]
				#print " "
				#print "Masses:" + str(zip(ID[sortList],mass_part_solar[sortList]))
				#print "Masses:" + str(mass_part_solar[sortList])

#		else:
#			continue
#	
	return numpy.array(times), numpy.array(masses), number_of_stars, start_time, end_time


def plot(dir_label):
	qtime, qmass, number_stars, start_time, end_time = run_quad(dir_label)
	print qtime, qmass, number_stars
	print 'last time in octant', qtime[-1]
	#times = numpy.arange(0, bins, 1) * (end_time-start_time)/bins + start_time
	#print times.size, start_time, end_time, end_time - start_time
	#import scipy.interpolate as si
	#totalMass = numpy.zeros( bins) 
	#totalnumber_stars = numpy.zeros(bins)
	boolMass = qmass > 0.1
	print boolMass
	print " "
	print qtime
	print qtime - start_time
	print end_time, start_time
	print qtime[0], start_time, qmass[0]

	pl.loglog( (qtime[boolMass] -start_time)/3.15e7, qmass[boolMass], label= r'N_J = ' + str(dir_label))
	#pl.loglog( (qtime -start_time)/3.15e7, qmass, label= r'N_J = ' + str(dir_label))
	#pl.loglog( (times[boolMass] -start_time)/3.15e7, qmass[boolMass], label= r'N_J = ' + str(dir_label))
	pl.loglog( numpy.arange(2, 10, 0.1)*1e5, 1.*(numpy.arange(2,10,0.1))**2, ls="--") 


# Constants #
mp = 1.6e-24
pi = numpy.pi
parsec = 3e18
sec_in_yr = numpy.pi* 1e7
Msun = 2e33
G = 6.67e-8
c_s = 2.64e4

units_time = 3.87201e3

bins = 1001
qtimes = []
qmasses = []
totnumber_stars = []
#quad_list = [2]
#for j in range(0, 1) : 
#	#quad = quad_list[j]
#	#print quad
#	qtime, qmass, number_stars = run_quad()
#	qtimes.append(qtime)
#	qmasses.append(qmass)
#	totnumber_stars.append(number_stars)
#	print 'time of 1st particle', qtime[0]
#	print 'last time in octant', qtime[-1]

home = "/home/m/murray/dwmurray/" 
os.chdir(home + "scratch/test-ramses/hd_32jeans_1024_v1")
plot('32')
##os.chdir(home + "test-ramses/jeans_32c")
##plot('32c')
#os.chdir(home + "test-ramses/jeans_16")
#plot('16')
#os.chdir(home + "test-ramses/jeans_8")
#plot('8')
#os.chdir(home + "test-ramses/jeans_4")
#plot('4')
#pl.ylim(1.2e0, 1e3)
pl.xlim(1.e4, 2e6)
pl.ylabel('$M_*$ $(M_\odot)$', fontsize = 25)
pl.xlabel( "$t-t_*$ $({\\rm yrs})$", fontsize = 25)

pl.legend(loc=0)
pl.savefig(home + "scratch/Total_mass_over_time.pdf")

#
#qtime, qmass, number_stars, start_time, end_time = run_quad()
#print qtime, qmass, number_stars
##qtimes.append(qtime)
##qmasses.append(qmass)
##totnumber_stars.append(number_stars)
##print 'time of 1st particle', qtime[0]
#print 'last time in octant', qtime[-1]
#
# find the start, end time
##start_time = qtime[0] #1e30
##end_time = qtime[-1]#1e30
##for qtime  in qtimes : 
##	if( qtime[0] < start_time) :
##		start_time = qtime[0]
##		print 'particle start time of the box', start_time
##	if( qtime[-1] < end_time) : 
##		end_time = qtime[-1]
##		print 'The earliest end time', end_time
#
#bins = 1001
#times = numpy.arange(0, bins, 1) * (end_time-start_time)/bins + start_time
##times = numpy.arange( start_time, end_time, (end_time - start_time)/bins)
#print times.size, start_time, end_time, end_time - start_time
#import scipy.interpolate as si
#totalMass = numpy.zeros( bins) 
#totalnumber_stars = numpy.zeros(bins)
##for qtime, qmass, number_stars  in zip( qtimes, qmasses, totnumber_stars) :
##	
##	func = si.interp1d(qtime, qmass, bounds_error = False, fill_value = 0e0, kind=3) 
##	totalMass = totalMass +  func(times)
##	func2 = si.interp1d(qtime, number_stars, bounds_error = False, fill_value = 0e0, kind=3) 
##	totalnumber_stars = totalnumber_stars  + func2(times)
##
##for i in qmass:
#	#print i
#	#print totalMass[i]
#
#times -= start_time
#times_yr = times / (numpy.pi*1.e7)
#times_hkyr = times_yr /(1.e5)
#
#boolMass = qmass > 0.1
#print boolMass
##print len(times_yr), len(totalnumber_stars)
##print totalnumber_stars
##pl.loglog( times_yr, totalMass) 
#print " "
#print start_time
##print qtime
##print (qtime -start_time)/3.15e7
##pl.loglog( (qtime -start_time)/3.15e7, qmass)
#pl.loglog( (qtime[boolMass] -start_time)/3.15e7, qmass[boolMass])
#pl.loglog( numpy.arange(4, 15, 0.1)*1e8, 3.*(numpy.arange(4,15,0.1))**2, ls="--") 
##pl.loglog( (times[boolMass] -start_time)/3.15e7, qmass[boolMass])#, ls='.-') 
##pl.loglog( times_yr, totalnumber_stars) 
##pl.loglog( times_hkyr, (times_hkyr/1e12)**2) 
##pl.loglog( numpy.arange(2, 5, 0.1)*1e4, 4.*(numpy.arange(2,5,0.1))**1, ls="-.-") 
##pl.ylim(1.2e0, 3e2)
##pl.xlim(3.e4, 6e5)
##pl.ylim(1.2e-1, 3e2)
##pl.xlim(1.e3, 6e5)
##pl.xticks(fontsize=25)
##pl.yticks(fontsize=25)
##pl.tick_params('both',length=8, width=1, which='minor')
##pl.tick_params('both',length=10, width=1.5, which='major')
##pl.gcf().subplots_adjust(bottom=0.15)
#pl.ylabel('$M_*$ $(M_\odot)$', fontsize = 25)
#pl.xlabel( "$t-t_*$ $({\\rm yrs})$", fontsize = 25)
#pl.savefig("/home/m/murray/dwmurray/scratch/Total_mass_over_time.pdf")

#for t,m in zip((qtime-start_time)/3.15e7, qmass) :
#	print "{0} {1}".format(t,m)
