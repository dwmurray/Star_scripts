import shutil
import argparse
import os
import numpy as np
#from yt.mods import * # depricated in yt3
import yt
import matplotlib.colorbar as cb
import glob
import sys
#from PIL import Image

def annotate_particles( yt_object, sinkfile, plot_type="density", width=(16, "pc")) :
	print sinkfile
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
		mass, rstar, xstar, ystar, zstar = np.loadtxt( sinkfile, usecols=[1,2,3,4,5], unpack=True, skiprows=3, comments="=")
		if( mass.size > 0 ) : 
			ax.scatter(xstar/pc+xlim[0], ystar/pc+ylim[0], color="black")
			# reset the limits
			ax.set_xlim(xlim)
			ax.set_ylim(ylim)
	except ValueError : 
		1+1
	os.remove(junkfile)


def particle_ID_locations(sinkfile, current_item) :
	if( not os.path.isfile( sinkfile)) :
		return
	try : 
		ID, mass, rstar, xstar, ystar, zstar = np.loadtxt( sinkfile, usecols=[0,1,2,3,4,5], unpack=True, skiprows=3, comments="=")
	except ValueError : 
		1+1
	print ID[current_item], mass[current_item], xstar[current_item], ystar[current_item], zstar[current_item]
	return ID[current_item], mass[current_item], xstar[current_item], ystar[current_item], zstar[current_item]

def Fullprojection_plotting(pf):
        plot_out_prefix = 'movieframe'
	print "Doing a projection plot of the Full Box."
	p = yt.ProjectionPlot(pf, "z", "density")
	fileout="{0}_{1:05d}_fullprojection.{2}".format(plot_out_prefix, i, out_format)
	if(withParticles) :
		print 'Annotating particles now'
		annotate_particles(p, sinkfile)
		fileout="{0}_{1:05d}_fullproject_particle.{2}".format(plot_out_prefix, i, out_format)
	p.set_font({'size':25})
	p.set_zlim("density", 1e-3 , 1e0)
	print fileout
	p.save(fileout)


def Projection_plotting(pf, ParticleID, xc, yc, zc, zoom_width=2.0):
        plot_out_prefix = 'movieframe'
	#plot_field = 'velocity' + plot_axis
	plot_field = 'density'
	fileout="{0}_{1:05d}_{2}_{3}_{4}pc_{5}.{6}".format(plot_out_prefix, i, plot_axis, plot_field, zoom_width, ParticleID, out_format)
	print plot_axis
	print plot_field
	#p = yt.ProjectionPlot(pf, plot_axis, plot_field, weight_field = 'density', center = (xc, yc, zc), width = (zoom_width, 'pc'))
	p = yt.ProjectionPlot(pf, 'z', 'density', center = (xc, yc, zc), width = (zoom_width, 'pc'))
	print 'have projected sucessfully'
	#p.set_zlim("Density", 1e-23 , 1e-14)
	#pid = os.getpid()
	#p.set_font({'size':25})
#	if( withParticles ) :
#		annotate_particles(p, sinkfile)
	print fileout
	#The issue right now is in saving it
	p.save(fileout)


def Slice(pf, xc, yc, zc, zoom_width=2.0) :
        plot_out_prefix = 'movieframe'
	sp = pf.h.sphere([xc, yc, zc], (zoom_width + 0.5, "pc"))
	# Get the angular momentum vector for the sphere.
	L = sp.quantities.angular_momentum_vector()
	# Create an OffAxisSlicePlot on the object with the L vector as its normal
	plot_field = 'density'

	# What axis do we plot along?
	Axis_to_plot = L
#	p = yt.OffAxisSlicePlot(pf, Axis_to_plot, plot_field, center=sp.center, width=(zoom_width, "pc"))
	p = yt.OffAxisSlicePlot(pf, Axis_to_plot, plot_field, center=(xc, yc, zc), width=(zoom_width, "pc"))
	p.set_zlim("density", 1e-23,1e-14)
	if(withArrows):
		#p.annotate_velocity(factor=16)
		p.annotate_contour("density")
	pid = os.getpid()
	p.set_font({'size':25})
        fileout="{0}_{1:05d}_{2}_{3}pc_{4}.{5}".format(plot_out_prefix, i, plot_field, zoom_width, this_particle_id, out_format)
	print fileout
        p.save(name=fileout)


parser = argparse.ArgumentParser(description = "start number to end number, step size, zoom width, plot axis and optional particle ID")
parser.add_argument('start', metavar='N1', type=int)
parser.add_argument('end', metavar='N2', type=int)
parser.add_argument('step', metavar='N3', type=int)
parser.add_argument('zoom_width', metavar='N4', type=float)
parser.add_argument('plot_axis', metavar='axis_', type=str, nargs='?', default='z', help='Axis along which to project/slice')
parser.add_argument('ParticleID', metavar='N5', type=int, nargs='?', default=42, help='Particle ID you want to reduce.')
parser.add_argument('--particle', action='store_true')
parser.add_argument('--noparticle', action='store_true')
parser.add_argument('--location', action='store_true')
# What kind of visualization do we want?
parser.add_argument('--project', action='store_true')
parser.add_argument('--projectfull', action='store_true')
parser.add_argument('--slice', action='store_true')
parser.add_argument('--arrow', action='store_true')
parser.add_argument('--allparticles', action='store_true')
parser.add_argument('--pdf', action='store_true')
#parser.add_argument('--parsec', metavar='N4', action='store_true')
args = parser.parse_args()
zoom_width = args.zoom_width
withParticles=args.particle
withNoParticles = args.noparticle
withProjection=args.project
withProjectionFull=args.projectfull
withLocationSpecified=args.location
withParticleIDValue = args.ParticleID
withAllParticles = args.allparticles
withSlice=args.slice
withArrows=args.arrow
plot_axis = args.plot_axis
withPDF = args.pdf

prefix = 'output'

if (withPDF):
	out_format = 'pdf'
else:
	out_format = 'png'

for i in range(args.start,args.end,args.step) :
	file_plt_exist = glob.glob("{0}_{1:05d}".format(prefix, i))
	if not file_plt_exist:
		print 'File: "', "{0}_{1:05d}".format(prefix, i), '" does not exist, moving to the next file'
		continue

	print 'On File: {0}_{1:05d}'.format(prefix, i)
	file = '{0}_{1:05d}/info_{1:05d}.txt'.format(prefix, i)
	print file

        sinkfile = None
        if( args.particle) :
                sinkfile = "{0}_{1:05d}/sink_{1:05d}.info".format(prefix, i)

	print sinkfile
 	pf = yt.load(file)

	print type(sinkfile)
	if( not os.path.isfile( sinkfile)) :
		print 'sinkfile doesn" exist here'
		sys.exit()
	try : 
		ID = np.loadtxt( sinkfile, usecols=[0], unpack=True, skiprows=3, comments="=")
	except ValueError : 
		1+1
	
	if (withProjectionFull):
		Fullprojection_plotting(pf)
	elif (withProjection):
		for current_item in range(len(ID)):
			ParticleID, ParticleMass, xc, yc, zc = particle_ID_locations(sinkfile, current_item)
		#Proj_plotting(pf, zoom_width)
			Projection_plotting(pf, ParticleID, xc, yc, zc, zoom_width=2.0)
	elif (withSlice):
		for current_item in range(len(ID)):
#			ParticleID, Particlemass, ParticleX, ParticleY, ParticleZ = particle_ID_locations(sinkfile, current_item)
			Particle_ID, ParticleMass, xc, yc, zc = particle_ID_locations(sinkfile, current_item)
			if Particle_ID == withParticleIDValue:
				Slice(pf, xc, yc, zc, zoom_width)
		#slice( file, xc, yc, zc, zoom_width)
	else:
		print 'need to specify what plotting, slice or projection.'
			
print "Finished, closing up shop"


	#project( file, moviefile, withParticles=args.particle)
	#mslice( file, moviefile)
        #slice( file, moviefile)
	#slice( file, moviefile,field="accx")
	#multislice( file, field="Density", fileout=moviefile)
	#multislice( "one_proc.hdf5", field="accx", fileout=moviefile)
	#ray_trace(file, moviefile)
	#histogram( file, field="Density", fileout=moviefile)
