import shutil
import argparse
import os
import numpy as np
from yt.mods import *
import matplotlib.colorbar as cb
from PIL import Image

import glob
import sys
import math
import h5py


#def project( filein, xc, yc, zc, zoom_width, fileout = 'junk.png', withParticles=False) :
#	if(withProjectionFull):
#		print "Doing a projection plot of the Full Box."
#		p = ProjectionPlot(pf, "z", "Density")
#	if not (withProjectionFull):
#		p = ProjectionPlot(pf, "z", "Density", center = (xc, yc, zc), width = (zoom_width, 'pc'))
#        p.set_zlim("Density", 1e-3,1e0)
#	if(withParticles) :
#		p.annotate_particles(1.0/pf['unitary'], p_size=10.0)
#        pid = os.getpid()
#        p.save("junk{0:06d}".format(pid))
#        shutil.move("junk{0:06d}_Projection_z_Density.png".format(pid),fileout)
#
def Proj_plotting(pf, xc, yc, zc, zoom_width, quad, i, this_particle_id):
        plot_out_prefix = 'movieframe'
	if(withProjectionFull):
		print "Doing a projection plot of the Full Box."
		p = ProjectionPlot(pf, "z", "Density")
		fileout="{0}_{1}_{2:04d}_fullprojection.{3}".format(plot_out_prefix, quad, i, out_format)
		if(withParticles) :
			p.annotate_particles(1.0/pf['unitary'], p_size=10.0)
			fileout="{0}_{1}_{2:04d}_fullproject_particle.{3}".format(plot_out_prefix, quad, i, out_format)
		p.set_font({'size':25})
		p.set_zlim("Density", 1e-3 , 1e0)
		p.set_colorbar_label("$\\Sigma\\,({\\rm g\\,cm}^{-2})$",fontsize=25)
		#plot = p.plots[field]
		#plot.axes = grid[i].axes
		#p.axes = 'Surface Density'
		#pid = os.getpid()
		print fileout
		p.save(fileout)
		#p.save("junk{0:06d}".format(pid))
		#shutil.move("junk{0:06d}_Projection_{1}_Density.png".format(pid, plot_axis),fileout)
	if not (withProjectionFull):
		#plot_field = plot_axis + '-velocity'
		plot_field = 'Density'
		fileout="{0}_{1}_{2:04d}_{3}_{4}_{5}pc_{6:03d}.{7}".format(plot_out_prefix, quad, i, plot_axis, plot_field, zoom_width, this_particle_id, out_format)
		p_plot = ProjectionPlot(pf, plot_axis, plot_field, center = (xc, yc, zc), width = (zoom_width, 'pc'))
		pid = os.getpid()
		p_plot.set_font({'size':25})
		#p_plot.set_zlim("Density", 1e-23 , 1e-14)
		if( withParticles ) :
			p_plot.annotate_particles(1.0/pf['unitary'], p_size=10.0)
		print fileout
		p_plot.set_colorbar_label("$\\Sigma\\,({\\rm g\\,cm}^{-2})$",fontsize=25)
		p_plot.save(fileout)
#		p_plot.save("junk{0:06d}".format(pid))
#		shutil.move("junk{0:06d}_Projection_{1}_{2}_Density.png".format(pid, plot_axis, plot_field),fileout)


def slice(filein, xc=0.0, yc=0.0, zc=0.0, zoom_width = 2.0, fileout='junk.png', withParticles=False, field="Density") :
	pf = load(file)
        plot_out_prefix = 'movieframe'
	#sp = pf.h.sphere([xc, yc, zc], (zoom_width + 0.5, "pc"))
	sp = pf.h.sphere([xc, yc, zc], (0.01, "pc"))
	# Get the angular momentum vector for the sphere.
	L = sp.quantities["AngularMomentumVector"]()
	# Find vectors orthogonal to L to plot side on view of disk.
	orthog_to_L_1 = np.array([L[2], L[2], -L[0]-L[1]])
	orthog_to_L_2 = np.array([L[1], -L[0]-L[2], L[1]])
	print np.dot(L, orthog_to_L_1)
	print np.dot(L, orthog_to_L_2)
	# What axis do we plot along?
	Axis_to_plot = L
	#Axis_to_plot = orthog_to_L_1
	#Axis_to_plot = orthog_to_L_2

	# Create an OffAxisSlicePlot on the object with the L vector as its normal
	plot_field = 'Density'
	#plot_field = 'VelocityMagnitude'

 	p = OffAxisSlicePlot(pf, Axis_to_plot, plot_field, sp.center, (zoom_width, "pc"))
	p.set_zlim("Density", 1e-23,1e-14)
	#p.set_zlim("Density", 1e-20,1e-17)
	if(withArrows):
		#v_vector_bulk=pf.h.disk([xc, yc, zc], Axis_to_plot, (1e-2, "pc"), (0.001, "pc")).quantities["BulkVelocity"]()
		v_vector_bulk=pf.h.sphere([xc, yc, zc], (2e-2, "pc")).quantities["BulkVelocity"]()
		p = OffAxisSlicePlot(pf, Axis_to_plot, plot_field, sp.center, (zoom_width, "pc"), field_parameters={"bulk_velocity": v_vector_bulk})
		p.set_zlim("Density", 1e-23,1e-14)
		p.annotate_cquiver('CuttingPlaneVelocityX', 'CuttingPlaneVelocityY', 12)
		#p.annotate_contour("Density")
	pid = os.getpid()
	#p.set_font({'family':'sans-serif', 'style':'italic', 'weight':'bold', 'size':24, 'color':'blue'})
	p.set_font({'size':25})
	#p.set_colorbar_label("$\\Sigma\\,({\\rm g\\,cm}^{-2})$",fontsize=25)
        fileout="{0}_{1}_{2:04d}_{3}_{4}pc_{5}.{6}".format(plot_out_prefix, quad, i, plot_field, zoom_width, this_particle_id, out_format)
	print fileout
        p.save(name=fileout)#, suffix='.pdf')

def mslice( filein, fileout = 'junk.png',field="Density") :
	newImage = Image.new("RGB", (1070,2000))
	pf = load(file)
	print moviefile
	center = na.array([0.0, 0.0, -5.795e16])
#	center = na.array([0.0, 0.0, 0e12])
	#p = SlicePlot(pf, "z", "dens", center=center )
	p = SlicePlot(pf, "x", field, center="max", width=(4,"pc") )
	p.set_zlim("Density", 1e-23,1e-19)
	p.annotate_velocity(factor=16)
	
	pid = os.getpid()
        p.save("junk{0:06d}".format(pid))

	p = SlicePlot(pf, "y", field, center="max", width=(4,"pc") )
	p.set_zlim("Density", 1e-23,1e-19)
	p.annotate_velocity(factor=16)
        p.save("junk{0:06d}".format(pid))

	file1 = "junk{0:06d}_Slice_x_{1}.png".format(pid,field)
	file2 = "junk{0:06d}_Slice_y_{1}.png".format(pid,field)
	image1 = Image.open(file1)
	image2 = Image.open(file2)

	newImage.paste( image1, (0,0))
	newImage.paste( image2, (0,1000))
	newImage.save(fileout)

	os.remove(file1)
	os.remove(file2)

def multislice( filein, field="Density", fileout='junk.png') :
	pf = load(filein)
	orient = "horizontal"
	
	fig, axes, colorbars = get_multi_plot(2, 1, colorbar=orient, bw=4)
	#pc = PlotCollection(pf, [0., 0., 0.])
	pc = PlotCollection(pf, [2.4e19,2.4e19,2.4e19])

	p = pc.add_slice( field, "x", axes= axes[0][0], use_colorbar=False)
	p.set_cmap("bds_highcontrast")
	#p.modify["velocity"](factor=16)
	p.annotate_velocity(factor=16)

	p = pc.add_slice( field, "y", axes= axes[0][1], use_colorbar=False)
	p.set_cmap("bds_highcontrast")
	#p.show_velocity(factor=16)
	#p = pc.add_slice( field, "z", axes= axes[0][2], use_colorbar=False)
	#p.set_cmap("bds_highcontrast")
	#p.show_velocity(factor=16)

	for p, cax in zip(pc.plots, colorbars) : 
		cbar = cb.Colorbar(cax, p.image, orientation=orient)
		p.colorbar = cbar
		p._autoset_label()

	fig.savefig(fileout)
	

def ray_trace( filein, fileout='junk.png') : 
	pf = load( filein)
	pf.h
	field = "Density"
	pf.field_info[field].take_log = True
	
	dd = pf.h.all_data()
	mi, ma = dd.quantities["Extrema"](field)[0]

	mi, ma = na.log10(mi), na.log10(ma)

	tf = ColorTransferFunction((mi, ma))

	c = pf.domain_center
	L = na.array([0.0, 1.0, 0.0])
	W = 1.0/pf["unitary"]
	N = 512

	cam = pf.h.camera(c, L, W, N, tf)
	tf.add_layers(10, 0.01, colormap = "RdBu_r")
	cam.snapshot(fileout)

def histogram( filein, field="Density", fileout='junk.png', xbounds=[1e-24, 1e-18]) :
        pf = load(filein)
	pc = PlotCollection(pf)
	data = pf.h.all_data()
	if( xbounds == None) : 
		pc.add_profile_object(data, [field, "CellMassMsun"], weight=None)
	else :
		pc.add_profile_object(data, [field, "CellMassMsun"], weight=None, x_bounds=xbounds)
	
	pid = os.getpid()
        pc.save("junk{0:06d}".format(pid))
	shutil.move("junk{0:06d}_Profile1D_0_{1}_CellMassMsun.png".format(pid,field),fileout)



def find_particle_position( filein):
	pf = load("{0}{1:04d}".format( plt_prefix, i))

	dd = pf.h.all_data()
	xp = dd["particle_posx"]
	yp = dd["particle_posy"]
	zp = dd["particle_posz"]
	partMass = dd["ParticleMassMsun"]


parser = argparse.ArgumentParser(description = "start number to end number")

parser.add_argument('start', metavar='N1', type=int)
parser.add_argument('end', metavar='N2', type=int)
parser.add_argument('step', metavar='N3', type=int)
parser.add_argument('zoom_width', metavar='N4', type=float)
parser.add_argument('plot_axis', metavar='axis_', type=str, nargs='?', default='z', help='Axis along which to project/slice')
parser.add_argument('ParticleID', metavar='N5', type=int, nargs='?', default=42, help='Particle ID you want to reduce.')
parser.add_argument('--particle', action='store_true')
parser.add_argument('--chk', action='store_true')
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
withCheckpoint = args.chk
withPDF = args.pdf

if (withCheckpoint) :
	prefix = "BB_hdf5_chk_"	
else:
	prefix = "BB_hdf5_plt_cnt_"
quad = os.getcwd()[-5:]

#plt_prefix = "BB_hdf5_plt_cnt_"
#part_prefix = "BB_hdf5_part_"
if (withPDF):
	out_format = 'pdf'
else:
	out_format = 'png'

if zoom_width > 8.0:
	print "The width of the box extends beyond the data values"
	print "that this quadrant calculated with AMR."
	print "System now exiting."
	sys.exit()

for i in range(args.start,args.end,args.step) :
	fn_plt = prefix+"{0:04d}".format(i)
	#fn_part = part_prefix+"{0:04d}".format(i)
	file_plt_exist = glob.glob(fn_plt)
	if not file_plt_exist:
		print 'File: "', fn_plt, '" does not exist, moving to the next file'
		continue
	#file_part_exist = glob.glob(fn_part)
	#if not file_part_exist:
	#	print 'File: "', fn_part, '" does not exist, moving to the next file'
	#	continue
	print 'On File: {0}_{1:04d}'.format(prefix, i)

 	pf = load("{0}{1:04d}".format(prefix, i))
	dd = pf.h.all_data()
	if (withNoParticles):
		max= dd.quantities["MaxLocation"]("Density")
		maxDens = max[0]/3e-22
		maxLoc = numpy.array(max[2:5])/3e18
		xc = max[2]
		yc = max[3]
		zc = max[4]
		#xc = 3.01906738281e+19
		#yc = 2.04562988281e+19
		#zc = 1.02202148438e+19
		#print xc
		this_particle_id = withParticleIDValue
		if (withLocationSpecified):
			this_particle_id = withParticleIDValue
			this_particle_id = str(this_particle_id) + '.0'
			if len(str(i)) == 2:
				four_digit_file_number = '00' + str(i)
			elif len(str(i)) == 3:
				four_digit_file_number = '0' + str(i)
			elif len(str(i)) == 4:
				four_digit_file_number = str(i)

			loaded_four_digit_file_number, loaded_withparticleidvalue, xc_search, yc_search, zc_search, this_creation_time, old_pulled_particle_mass, old_pulled_current_time = np.loadtxt("particle_{0}_location.txt".format(withParticleIDValue), unpack=True)

			#print loaded_four_digit_file_number
			#print loaded_withparticleidvalue
			for value in range(len(loaded_four_digit_file_number)):
				pulled_four_digit_file_number = loaded_four_digit_file_number[value]
				pulled_withparticleidvalue = loaded_withparticleidvalue[value]
				pulled_four_digit_file_number = str(pulled_four_digit_file_number)[:-2]
				if len(str(pulled_four_digit_file_number)) == 2:
					pulled_four_digit_file_number = '00' + str(pulled_four_digit_file_number)
				elif len(str(pulled_four_digit_file_number)) == 3:
					pulled_four_digit_file_number = '0' + str(pulled_four_digit_file_number)
				elif len(str(pulled_four_digit_file_number)) == 4:
					pulled_four_digit_file_number = str(pulled_four_digit_file_number)
				else:
					print 'something is wrong with the file number'
					sys.exit()

				print four_digit_file_number, pulled_four_digit_file_number
				print pulled_withparticleidvalue, withParticleIDValue
				if four_digit_file_number == pulled_four_digit_file_number:
					xc = xc_search[value]
					yc = yc_search[value]
					zc = zc_search[value]
					print xc, yc, zc
					break
				else:
					print 'The file numbers do not match'
					continue
		file = prefix+"{0:04d}".format(i)
		print file
		if (withProjection) or (withProjectionFull):
			moviefile = "movie_{0}_frame_{1:04d}_prjct_{2}pc_{3}.{4}".format(quad, i, zoom_width, this_particle_id, out_format)
			#project( file, xc, yc, zc, zoom_width, moviefile, withParticles)
			Proj_plotting(pf, xc, yc, zc, zoom_width, quad, i, this_particle_id)
		elif (withSlice):
			moviefile = "movie_{0}_frame_{1:04d}_slice_{2}pc_{3}.{4}".format(quad, i, zoom_width, this_particle_id, out_format)
			slice( file, xc, yc, zc, zoom_width, moviefile, withParticles)
		else:
			print 'need to specify what plotting, slice or projection.'
			
	else:
		xp = dd["particle_posx"]
		yp = dd["particle_posy"]
		zp = dd["particle_posz"]
		particle_ID = dd["particle_index"]
		#partMass = dd["ParticleMassMsun"]
		for j in range(xp.size) :
			print particle_ID[j], withParticleIDValue
			#print type(particle_ID[j]), type(withParticleIDValue)
			if particle_ID[j] == withParticleIDValue or (withAllParticles):
				xc = xp[j]
				yc = yp[j]
				zc = zp[j]
				this_particle_id = int(particle_ID[j])
				#print type(this_particle_id), this_particle_id
				file = prefix+"{0:04d}".format(i)
				print file
				print "Currently on particle ", j + 1, "of ", xp.size
				if (withProjection) or (withProjectionFull):
					moviefile = "movie_{0}_frame_{1:04d}_prjct_{2}pc_{3}.{4}".format(quad, i, zoom_width, this_particle_id, out_format)
			#project( file, xc, yc, zc, zoom_width, moviefile, withParticles)
					Proj_plotting(pf, xc, yc, zc, zoom_width, quad, i, this_particle_id)
				elif (withSlice):
					moviefile = "movie_{0}_frame_{1:04d}_slice_{2}pc_{3}.{4}".format(quad, i, zoom_width, this_particle_id, out_format)
					slice( file, xc, yc, zc, zoom_width, moviefile, withParticles)
				else:
					print 'did not specify a plot to do, slice or projection'

print "Finished, closing up shop"
	#project( file, moviefile, withParticles=args.particle)
	#mslice( file, moviefile)
        #slice( file, moviefile)
	#slice( file, moviefile,field="accx")
	#multislice( file, field="Density", fileout=moviefile)
	#multislice( "one_proc.hdf5", field="accx", fileout=moviefile)
	#ray_trace(file, moviefile)
	#histogram( file, field="Density", fileout=moviefile)
