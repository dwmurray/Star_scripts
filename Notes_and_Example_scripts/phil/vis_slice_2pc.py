import shutil
import argparse
import os
from yt.mods import *
import matplotlib.colorbar as cb
from PIL import Image

def project( filein, fileout = 'junk.png', withParticles=False) :
	pf = load(file)
	print moviefile
	p = ProjectionPlot(pf, "z", "Density")
	p.set_zlim("Density", 1e-3,1e0)
	if( withParticles ) : 
		p.annotate_particles(1.0/pf['unitary'], p_size=10.0)
	pid = os.getpid()
        p.save("junk{0:06d}".format(pid))
	shutil.move("junk{0:06d}_Projection_z_Density.png".format(pid),fileout)

def slice( filein, fileout = 'junk.png',field="Density") :
	pf = load(file)
	print moviefile
	print xc
	print yc
	print zc

	#sp = pf.h.sphere("max", (1e-2, "pc"))
	#sp = pf.h.sphere([xc, yc, zc], (1e-2, "pc"))
	sp = pf.h.sphere([xc, yc, zc], (1e0, "pc"))

	# Get the angular momentum vector for the sphere.
	L = sp.quantities["AngularMomentumVector"]()

	# Create an OffAxisSlicePlot on the object with the L vector as its normal
	p = OffAxisSlicePlot(pf, L, "Density", sp.center, (5e-1, "pc"))
        #p = OffAxisSlicePlot(pf, L, "Density", sp.center, (1e-1, "pc"))
	p.set_zlim("Density", 1e-23,1e-14)
	#p.annotate_velocity(factor=16)
	
	pid = os.getpid()
        p.save("junk{0:06d}".format(pid))
	shutil.move("junk{0:06d}_OffAxisSlice_{1}.png".format(pid,field),fileout)

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

	
prefix = "BB_hdf5_plt_cnt_"
#prefix = "BB_hdf5_chk_"

parser = argparse.ArgumentParser(description = "start number to end number")

parser.add_argument('start', metavar='N1', type=int)
parser.add_argument('end', metavar='N2', type=int)
parser.add_argument('step', metavar='N3', type=int)
parser.add_argument('--particle', action='store_true')

#parser.add_argument('--parsec', metavar='N4', action='store_true')

args = parser.parse_args()

quad = os.getcwd()[-5:]

for i in range(args.start,args.end,args.step) :

        pf = load("{0}{1:04d}".format(prefix, i))

        dd = pf.h.all_data()
        xp = dd["particle_posx"]
        yp = dd["particle_posy"]
        zp = dd["particle_posz"]
        partMass = dd["ParticleMassMsun"]

	for j in range(xp.size) :
                xc = xp[j]
                yc = yp[j]
                zc = zp[j]
#                getRadialProfile(pf,xc,yc,zc, fileout="{0}{1:04d}_part_{2:03d}.out".format( out_prefix, i, j), radiusSphere=3.0, particleMass =partMass[j])

		file = prefix+"{0:04d}".format(i)
		print file
		moviefile = "movie_{0}_frame_{1:04d}_part_{2:04d}_05pc.png".format(quad, i, j)
	#project( file, moviefile, withParticles=args.particle)
		slice( file, moviefile)
	#mslice( file, moviefile)

	#slice( file, moviefile,field="accx")
	#multislice( file, field="Density", fileout=moviefile)
	#multislice( "one_proc.hdf5", field="accx", fileout=moviefile)
	#ray_trace(file, moviefile)
	#histogram( file, field="Density", fileout=moviefile)
