import shutil
import matplotlib 
matplotlib.use("Agg")
import argparse
import os
import yt
import numpy
import matplotlib.pyplot as pl
import math

def annotate_particles( yt_object, sinkfile, plot_type="density", width=(16, "pc")) :
	if( not os.path.isfile( sinkfile)) :
		return
	# get axis 
	ax = yt_object.plots[plot_type].axes
	xlim = ax.get_xlim()
	ylim = ax.get_ylim()
	# save the plot to execute first
	pid = os.getpid()
	junkfile = "junk{0}.png".format(pid)
	yt_object.save(junkfile)
	pc = 3.0857e18
	mass, xstar, ystar, zstar = numpy.loadtxt( sinkfile, usecols=[1,2,3,4], unpack=True, skiprows=3, comments="=")
	if( mass.size > 0 ) : 
		ax.scatter(xstar/pc+xlim[0], ystar/pc+ylim[0], color="black")
		# reset the limits
		ax.set_xlim(xlim)
		ax.set_ylim(ylim)
	os.remove(junkfile)
	 	

def project( filein, fileout = 'junk.png', sinkfile=None) :
	pf = yt.load(file)
	print moviefile
	p = yt.ProjectionPlot(pf, "z", "density")
	p.set_zlim("density", 1e-3,1e0)
	if( not sinkfile == None) : 
		annotate_particles(p, sinkfile) 
		
        p.save(fileout)
	


parser = argparse.ArgumentParser(description = "start number to end number")

parser.add_argument('start', metavar='N1', type=int)
parser.add_argument('end', metavar='N2', type=int)
parser.add_argument('step', metavar='N3', type=int)
args = parser.parse_args()

#dirs = ["hires/jet", "hires/nojet", "medres/jet", "medres/nojet"]
dirs = ["highres/jet", "highres/nojet", "medres/jet", "medres/nojet"]
labels = ["$16K^3$ jet", "no jet", "$8K^3$ jet", "no jet"]
ltypes = ["solid", "dashed", "solid", "dashed"]
lweights = [4, 4, 2, 2]
for dir, label, ltype, lw in zip(dirs,labels,ltypes,lweights): 
	first_time = True
	t_start = 0.
	tarray = []
	mtot = []
	prefix = "output_"
	for i in range(args.start,args.end,args.step) :	
		file = "{0}/{1}".format(dir,prefix)+"{0:05d}/info_{0:05d}.txt".format(i)
		sinkfile = "{0}/{1}".format(dir,prefix)+"{0:05d}/sink_{0:05d}.info".format(i)
                if( not os.path.isfile(sinkfile)) : 
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
	#slice( file, moviefile)
	#mslice( file, moviefile)

	#slice( file, moviefile,field="accx")
	#multislice( file, field="Density", fileout=moviefile)
	#multislice( "one_proc.hdf5", field="accx", fileout=moviefile)
	#ray_trace(file, moviefile)
	#histogram( file, field="Density", fileout=moviefile)

	gasMass = 4.8e19**3 * 3e-22 / 2e33
	tarray = numpy.log10(tarray) 
	mtot = numpy.log10(mtot) - math.log10(gasMass)
	mtot = 1e1**mtot
	tarray = 1e1**tarray/3.1558e7
	for t, m in zip(tarray, mtot) :
		print t, m
	
	lm = numpy.log10(mtot)
        lt = numpy.log10(tarray)

        p = numpy.polyfit( lt[lt>5.5], lm[lt>5.5], 1)
	pl.loglog( tarray, mtot, linewidth=lw,label="{0} $\\alpha={1:2.1f}$".format(label, p[0]),ls=ltype)
        lt =  numpy.arange(5.5, 6.1, 0.1 )
	pl.loglog( 1e1**lt, 1e1**(p[0]*lt + p[1]), ls="dotted", linewidth=1,color="black")

pl.xlim(1e5, 2e6)
#arr =  numpy.arange(2e5, 3e6, 2e5)
#pl.loglog( arr, 0.1*(arr/1e6)**2, linewidth=2, ls="dotted",label="$\\propto (t-t_*)^2$")
pl.ylim(5e-4, 0.2)
pl.xlabel("$t-t_*$ [yrs]", fontsize=20)
pl.ylabel("$M_*/M_{\\rm tot}$", fontsize=20)
pl.legend(loc="best",fontsize=20)
pl.xticks(fontsize=17)
pl.yticks(fontsize=17)

pl.savefig("/home/m/murray/dwmurray/scratch/test.pdf")
