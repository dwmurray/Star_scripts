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


import h5py
import numpy
import os
import yt
import sys
import glob
import matplotlib.pyplot as pl

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
	print current_time / 3.15e7
	#current_time = 0.0
        return ID, mass, current_time

def simulation_info():
	times = []
	masses = []
	number_of_stars = []
	prefix = 'output_'
	start_time = 0
	calc_time = 0

	for fileIndex in range( 1, 500) : 
		file = prefix+"{0:05d}/info_{0:05d}.txt".format(fileIndex)
		sinkfile = prefix+"{0:05d}/sink_{0:05d}.info".format(fileIndex)
		print file, sinkfile
		if not glob.glob(file):
			continue
		if not glob.glob(sinkfile):
			continue
		print file, sinkfile
		ID, mass_part_solar, current_time = Obtain_particles(file, sinkfile, fileIndex)
		if ID == []:
			pass
		else :
			#end_time = current_time
			if start_time != 0:
				times.append(current_time)
				masses.append(mass_part_solar.sum())
				end_time = current_time
			else:
				start_time = current_time
				times.append(current_time)
				masses.append(mass_part_solar.sum())
				end_time = current_time
			try:
				number_of_stars.append(len(mass_part_solar))
			except TypeError:
				number_of_stars = [1]
			#sortList = numpy.argsort( mass_part_solar)[::-1]
			#print " "
			#print "Masses:" + str(zip(ID[sortList],mass_part_solar[sortList]))
#			return numpy.array(times), numpy.array(masses), number_of_stars, start_time, end_time
	return numpy.array(times), numpy.array(masses), number_of_stars, start_time, end_time


def plot():
	qtime, qmass, number_stars, start_time, end_time = simulation_info()
	#print qtime, qmass, number_stars
	#print 'last time in octant', qtime[-1]
	number_stars_at_mass = []
	mass_value =[]
	sortList = numpy.argsort(qmass)[::-1]
	qmass = qmass[sortList]
	print qmass
	for mass_limit in numpy.arange(1e-2, 1e2, 0.1):
		#print mass_value
		boolMass_above = qmass > mass_limit
		#print boolMass_above
		#print numpy.sum(boolMass_above)
		number_stars_at_mass.append(numpy.sum(boolMass_above))
		mass_value.append(mass_limit)
	
	withHistogram = True
	if withHistogram:
		log_mass = numpy.log10(qmass)
		bins = 20

		hist_tuple = numpy.histogram(log_mass,bins)
		Number_of_stars_per_bin = hist_tuple[0]+1.e-2 # to eliminate the log of zero issue if no starsin bin.
		x_bins = hist_tuple[1]
		max_numstars = numpy.amax(Number_of_stars_per_bin)
		test_mas = numpy.argmax(Number_of_stars_per_bin)
		print '# stars', Number_of_stars_per_bin
		print max_numstars
		print test_mas
		x_start_pt = x_bins[test_mas]
		print x_start_pt
		#sys.exit()

		boolMass = log_mass >= x_start_pt
		print boolMass
		#high_log_mass = numpy.log10(qmass[boolMass])
		high_log_mass = log_mass[boolMass]
		#print len(high_log_mass)
		#print log_mass

		Number_of_stars_per_bin = Number_of_stars_per_bin[test_mas:]
		x_bins = x_bins[test_mas:]

		#y_axis = numpy.histogram(high_log_mass,bins)
		#Number_of_stars_per_bin = y_axis[0]+1.e-2 # to eliminate the log of zero issue if no starsin bin.
		#x_bins = y_axis[1]
		x_mid = (x_bins[1:] + x_bins[:-1])/2.

		#Fit a straight line to the histogram above the midpoint of the curve
		#x_up = x_mid[bins/2:]
		y_up = Number_of_stars_per_bin#[bins/2:]
		print 'y_up = ', y_up
		x = x_mid
		y = numpy.log10(y_up)
		A = numpy.vstack([x, numpy.ones(len(x))]).T
		print len(y)
		print y
		print len(A)
		#print A
		m, c = numpy.linalg.lstsq(A,y)[0]
		print('slope, intercept=',m, c)

		y_fit = 10.**(m*x +c)

		pl.hist(log_mass, bins, histtype='step', color='black', log=True)
		#pl.plot(x_mid, Number_of_stars_per_bin, 'bo')
		#pl.plot(x_up, y_up, 'ro')
		pl.plot(x, y_fit, 'r')
		pl.xlabel(r'$\log_{10}(m)$', size=15)
		pl.ylabel(r'$\log N_{*}$', size=20)
#		pl.ylim(1.9e0, 2e1)
#		pl.xlim(-1.e0, 1.5e0)
		pl.savefig(home + "_imf_hist.pdf")
	else:
		pl.loglog(mass_value, number_stars_at_mass)
		#print boolMass
		#print " "
		#print qtime
		#print qtime - start_time
		#print end_time, start_time
		#print qtime[0], start_time, qmass[0]
		#pl.loglog( (qtime[boolMass] -start_time)/3.15e7, qmass[boolMass], label= r'N_J = ' + str(dir_label))
		#pl.loglog()
		#pl.loglog( numpy.arange(2, 10, 0.1)*1e5, 1.*(numpy.arange(2,10,0.1))**2, ls="--") 

		#pl.xlim(1.e-1, 1e2)
		#pl.ylabel('$N$ ', fontsize = 25)
		#pl.xlabel( "$M_*$ $(M_\odot)$", fontsize = 25)
		pl.savefig(home + "_imf.pdf")




# Constants #
mp = 1.6e-24
pi = numpy.pi
parsec = 3e18
sec_in_yr = numpy.pi* 1e7
Msun = 2e33
G = 6.67e-8
c_s = 2.64e4

units_time = 3.87201e3

qtimes = []
qmasses = []
totnumber_stars = []

#home = "/home/m/murray/dwmurray/" 
home = os.getcwd()#"/Users/dwmurray/Work/ramses-jet/"
print home
#os.chdir(home)
plot()


