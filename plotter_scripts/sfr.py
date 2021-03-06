import shutil
import matplotlib 
matplotlib.use("Agg")
import argparse
import os
import sys
import yt
import numpy
import matplotlib.pyplot as pl
import math
from matplotlib import rcParams

# Setup plot tick marks to match the other run of __ plots.
#pl.rc("text", usetex=True)
#pl.rc('font', family='serif')
#rcParams.update({'figure.autolayout': True})
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

parser = argparse.ArgumentParser(description = "start number to end number")

parser.add_argument('start', metavar='N1', type=int)
parser.add_argument('end', metavar='N2', type=int)
parser.add_argument('step', metavar='N3', type=int)
parser.add_argument('--kink', action='store_true')
parser.add_argument('--check', action='store_true')
args = parser.parse_args()

##dirs = ["hires/jet", "hires/nojet", "medres/jet", "medres/nojet"]
dirs = ["ultrares/jet", "ultrares/nojet", "highres/jet", "highres/nojet", "highres/jet_v0_1"]
#dirs = ["ultrares/testjet", "ultrares/testnojet", "highres/jet", "highres/nojet"]
labels = ["$32K^3$ 0.3 jet", "no jet", "$16K^3$ jet", "no jet", "$16K^3$ 0.1 jet"]
ltypes = ["solid", "dashed", "solid", "dashed", "solid"]
lweights = [4, 4, 2, 2, 2]

if args.check:
	dirs = ["highres/jet", "highres/nojet", "highres/jet_v0_1", "highres/jet_v0_5"]
	labels = ["$16K^3$ 0.3 jet (standard)", "no jet", "0.1 jet", "0.5 jet"]
	ltypes = ["solid", "dashed", "solid", "solid"]
	lweights = [2, 2, 2, 2]

for dir, label, ltype, lw in zip(dirs,labels,ltypes,lweights): 
	first_time = True
	t_start = 0.
	tarray = []
	mtot = []
	prefix = "output_"
	for i in range(args.start,args.end,args.step) :	
		file = "{0}/{1}".format(dir,prefix)+"{0:05d}/info_{0:05d}.txt".format(i)
		sinkfile = "{0}/{1}".format(dir,prefix)+"{0:05d}/sink_{0:05d}.info".format(i)
		print sinkfile
                if( not os.path.isfile(sinkfile)) : 
			continue
                if( not os.path.isfile(file)) : 
			continue
                print "opening {0}".format(sinkfile)
		mass = None
		try : 
			mass, xstar, ystar, zstar = numpy.loadtxt( sinkfile, usecols=[1,2,3,4], unpack=True, skiprows=3, comments="=")
			pf = yt.load(file)
			time = pf.current_time#.convert_to_cgs()
			if( first_time and mass.size > 0) :
				t_start = time
				first_time = False
			if( not first_time) :
				mtot.append(mass.sum())
				tarray.append(time-t_start)
                                #print "time = " + str(tarray[-1]/3.15e7)
		except ValueError : 
			print sinkfile + " empty" 
	gasMass = 4.8e19**3 * 3e-22 / 2e33
	tarray = numpy.log10(tarray) 
	mtot = numpy.log10(mtot) - math.log10(gasMass)
	mtot = 1e1**mtot
	tarray = 1e1**tarray/3.1558e7
	for t, m in zip(tarray, mtot) :
		print t, m
	
	lm = numpy.log10(mtot)
        lt = numpy.log10(tarray)

	lt_min = 5.5
	lt_max = 6.1
        p = numpy.polyfit( lt[lt>lt_min], lm[lt>lt_min], 1)

	if( args.kink):
		lt_kink = 5.9
		lt_mask = lt[lt>lt_min]
		lm_mask = lm[lt>lt_min]
		p = numpy.polyfit( lt_mask[lt_mask<lt_kink], lm_mask[lt_mask<lt_kink], 1)
		print p[deg]
		sys.exit()
		k = numpy.polyfit( lt[lt>=lt_kink], lm[lt>=lt_kink], 1)
		pl.loglog( tarray, mtot, linewidth=lw,label="{0} $\\alpha_1={1:2.2f}$ $\\alpha_2={2:2.2f}$".format(label, p[0], k[0]),ls=ltype)
		ltk =  numpy.arange(lt_kink, 6.2, 0.1 )
		pl.loglog( 1e1**ltk, 1e1**(k[0]*ltk + k[1]), ls="dotted", linewidth=1,color="black")
		lt =  numpy.arange(lt_min, lt_kink, 0.1 )
		pl.loglog( 1e1**lt, 1e1**(p[0]*lt + p[1]), ls="dotted", linewidth=1,color="black")
	else:
		pl.loglog( tarray, mtot, linewidth=lw,label="{0} $\\alpha_1={1:2.3f}$".format(label, p[0]),ls=ltype)
		#This is for plotting the fit, p. It will always be plotted.
		lt =  numpy.arange(lt_min, lt_max, 0.1 )
		pl.loglog( 1e1**lt, 1e1**(p[0]*lt + p[1]), ls="dotted", linewidth=1,color="black")


pl.xlim(5e4, 1e6)
pl.ylim(5e-5, 0.2)
#arr =  numpy.arange(2e5, 3e6, 2e5)
#pl.loglog( arr, 0.1*(arr/1e6)**2, linewidth=2, ls="dotted",label="$\\propto (t-t_*)^2$")

pl.xlabel("$t-t_*$ [yrs]", fontsize=20)
pl.ylabel("$M_*/M_{\\rm tot}$", fontsize=20)
pl.legend(loc="best",fontsize=15)
pl.xticks(fontsize=17)
pl.yticks(fontsize=17)

pl.savefig("/home/m/murray/dwmurray/scratch/ramses-jet/check_sfr.pdf")
#pl.savefig("sfr.pdf")
