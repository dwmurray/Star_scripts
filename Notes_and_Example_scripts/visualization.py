import shutil
import argparse
import os
import yt
import numpy

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
	try : 
		mass, xstar, ystar, zstar = numpy.loadtxt( sinkfile, usecols=[1,2,3,4], unpack=True, skiprows=3, comments="=")
		if( mass.size > 0 ) : 
			ax.scatter(xstar/pc+xlim[0], ystar/pc+ylim[0], color="black")
			# reset the limits
			ax.set_xlim(xlim)
			ax.set_ylim(ylim)
	except ValueError : 
		1+1
	os.remove(junkfile)
	 	

def project( filein, fileout = 'junk.png', sinkfile=None) :
	pf = yt.load(file)
	print moviefile
	p = yt.ProjectionPlot(pf, "z", "density")
	p.set_zlim("density", 1e-3,1e0)
	if( not sinkfile == None) : 
		annotate_particles(p, sinkfile) 
		
        p.save(fileout)
	
prefix = "output_"

parser = argparse.ArgumentParser(description = "start number to end number")

parser.add_argument('start', metavar='N1', type=int)
parser.add_argument('end', metavar='N2', type=int)
parser.add_argument('step', metavar='N3', type=int)
parser.add_argument('--particle', action='store_true')
args = parser.parse_args()

for i in range(args.start,args.end,args.step) :	
	file = prefix+"{0:05d}/info_{0:05d}.txt".format(i)
	print file
	moviefile = "movie_frame{0:04d}.png".format(i)
	sinkfile = None
	if( args.particle) : 
		sinkfile = prefix+"{0:05d}/sink_{0:05d}.info".format(i)
	project( file, moviefile, sinkfile)
	#slice( file, moviefile)
	#mslice( file, moviefile)

	#slice( file, moviefile,field="accx")
	#multislice( file, field="Density", fileout=moviefile)
	#multislice( "one_proc.hdf5", field="accx", fileout=moviefile)
	#ray_trace(file, moviefile)
	#histogram( file, field="Density", fileout=moviefile)
