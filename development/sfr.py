"""
This program finds the total mass in star particles over time in a simulation. It also finds the total number of stars over time.
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

################################################
################################################
################################################
################################################
################################################


parser = argparse.ArgumentParser(description = "start number to end number")

parser.add_argument('start', metavar='N1', type=int)
parser.add_argument('end', metavar='N2', type=int)
parser.add_argument('step', metavar='N3', type=int)
args = parser.parse_args()

# Constants #
mp = 1.6e-24
pi = np.pi
parsec = 3e18
#sec_in_yr = pi* 1e7
yr = 3.1558e7
Msun = 2e33
G = 6.67e-8
c_s = 2.64e4

global Init_matplotParams
Init_matplotParams = False


#Depending on if we are reducing the data on scinet or another location, choose the output file path accordingly
cwd = os.getcwd()
#output_location = cwd + '/python_output'

# Simulation Parameters #
box_length = 4.8e19
box_density = 3e-22
# Parameters #
bins = 1001
dirs = ["highres/jet", "highres/nojet", "medres/jet", "medres/nojet"]
labels = ["$16K^3$ jet", "no jet", "$8K^3$ jet", "no jet"]
ltypes = ["solid", "dashed", "solid", "dashed"]
lweights = [4, 4, 2, 2]
prefix = "output_"

Setup_plot_window()

for dir, label, ltype, lw in zip(dirs,labels,ltypes,lweights): 
	first_time = True
	t_start = 0.
	tarray = []
	mtot = []
	N_stars = []
	for i in range(args.start,args.end,args.step) :	
		file = "{0}/{1}".format(dir,prefix)+"{0:05d}/info_{0:05d}.txt".format(i)
		sinkfile = "{0}/{1}".format(dir,prefix)+"{0:05d}/sink_{0:05d}.info".format(i)
                if( not os.path.isfile(sinkfile)) : 
			continue
                print "opening {0}".format(sinkfile)
		mass = None
		try : 
			mass, rstar, xstar, ystar, zstar = np.loadtxt( sinkfile, usecols=[1,2,3,4,5], unpack=True, skiprows=3, comments="=")
			pf = yt.load(file)
			time = pf.current_time
			if( first_time and mass.size > 0) :
				t_start = time
				first_time = False
			if( not first_time) :
				mtot.append(mass.sum())
				tarray.append(time-t_start)
				N_stars.append(len(mass))

		except ValueError : 
			print sinkfile + " empty" 

	gasMass = box_length**3 * box_density / Msun
	tarray = np.log10(tarray) 
	mtot = np.log10(mtot) - math.log10(gasMass)
	mtot = 1e1**mtot
	tarray = 1e1**tarray/3.1558e7
	for t, m in zip(tarray, mtot) :
		print t, m
	
	lm = np.log10(mtot)
        lt = np.log10(tarray)

        p = np.polyfit( lt[lt>5.5], lm[lt>5.5], 1)
	plt.loglog( tarray, mtot, linewidth=lw,label="{0} $\\alpha={1:2.1f}$".format(label, p[0]),ls=ltype)
        lt =  np.arange(5.5, 6.1, 0.1 )
	plt.loglog( 1e1**lt, 1e1**(p[0]*lt + p[1]), ls="dotted", linewidth=1,color="black")


#arr =  np.arange(2e5, 3e6, 2e5)
#plt.loglog( arr, 0.1*(arr/1e6)**2, linewidth=2, ls="dotted",label="$\\propto (t-t_*)^2$")
plt.xlim(1e5, 2e6)
plt.ylim(5e-4, 0.2)
plt.xlabel("$t-t_*$ [yrs]", fontsize=20)
plt.ylabel("$M_*/M_{\\rm tot}$", fontsize=20)
plt.legend(loc="best",fontsize=20)
#plt.xticks(fontsize=17)
#plt.yticks(fontsize=17)

plt.savefig(cwd + "test.pdf")
