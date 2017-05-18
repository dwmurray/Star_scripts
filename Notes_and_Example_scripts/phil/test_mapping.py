import matplotlib
matplotlib.use("Agg")

from yt.mods import *
import matplotlib.pyplot as plt
import shutil
import math
import numpy
import random
import time
import datetime
import os
#from scipy import stats
import matplotlib.mlab as mlab
from scipy.stats import norm as scinorm

bins = 70
meanDens = 3e-22
i = 0
Msun = 1.99e33
G = 6.67e-8
parsec = 3.09e18
logRhoMin = -3.0

def getRadialProfile_yt(pf, xc, yc, zc, vxc, vyc, vzc, fileout="rad_profile.out", radiusMin=1e-3, radiusSphere=3.0, particleMass = 0.0) : 
	if (Sphere_Bulk):
		sp = pf.h.sphere((xc, yc, zc), radiusSphere/pf['pc'])
		sp2 = pf.h.sphere((xc, yc, zc), Bulk_sphere_radius/pf['pc'])
		print Bulk_sphere_radius
		bulk = sp2.quantities["BulkVelocity"]()
	if (Particle_Bulk):
		sp = pf.h.sphere((xc, yc, zc), radiusSphere/pf['pc'])
		#sp2 = pf.h.sphere((xc, yc, zc), 0.01/pf['pc'])
		bulk = numpy.array([vxc, vyc, vzc])
	if (NoParticle_Bulk):
		sp = pf.h.sphere((xc, yc, zc), 4.0/pf['pc'])
		bulk = sp.quantities["BulkVelocity"]()
	sp.set_field_parameter("bulk_velocity", bulk)
	rp = BinnedProfile1D( sp, bins, "Radiuspc", 1e1**logRhoMin, radiusSphere, log_space=True)
	rp.add_fields("Density", weight="CellVolume")
	maxDens = rp["Density"][0]
	
	rp.add_fields("RadialVelocity")  
	rp.add_fields("VelocityMagnitude")
	
	rp.add_fields("CellMassMsun", accumulation=True,weight=None)

	rbin = rp["Radiuspc"]
	vrbin = rp["RadialVelocity"]
	mTbin = rp["CellMassMsun"] + particleMass
	vKbin = numpy.sqrt(G*mTbin*Msun/(parsec*rbin))
	rhobin = rp["Density"]
	mdotbin = 4.0*3.141*rbin**2*vrbin*rhobin
	rp.add_fields("TangentialVelocity")
	vrmsbin = rp["TangentialVelocity"]
	rp.add_fields("TangentialVelocity", weight="CellVolume")
	vrmsnbin = rp["TangentialVelocity"]
	rp.add_fields("VelocityMagnitude")
	vmagbin = rp["VelocityMagnitude"]	
	rp.add_fields("VelocityMagnitude", weight="CellVolume")
	vmagnbin = rp["VelocityMagnitude"]

	rp.add_fields(["AngularMomentumX", "AngularMomentumY", "AngularMomentumZ"])

	norm = numpy.sqrt(rp["AngularMomentumX"]**2 + rp["AngularMomentumY"]**2 + rp["AngularMomentumZ"]**2)
	angXbin = rp["AngularMomentumX"]/norm
	angYbin = rp["AngularMomentumY"]/norm
	angZbin = rp["AngularMomentumZ"]/norm
	
	numpy.savetxt(fileout, zip(rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, rhobin, mdotbin, norm, angXbin, angYbin, angZbin), fmt="%15.9E")

def getRadialProfile_py(pf, xc, yc, zc, fileout="rad_profile.out", radiusMin=1e-3, radiusSphere=3.0, particleMass = 0.0) : 
	sp = pf.h.sphere((xc, yc, zc), radiusSphere/pf['pc'])

	x = sp["x"] - xc
	y = sp["y"] - yc
	z = sp["z"] - zc
	r = numpy.sqrt(x*x + y*y + z*z)
	# grams
	cellMass = sp["CellMassMsun"] * Msun
	# cgs
	dens = sp["Density"]
	# cm**3
	cellVolume = sp["CellVolume"]

	# These are in cm/s
	vx = sp["x-velocity"]
	vy = sp["y-velocity"]
	vz = sp["z-velocity"]

	# This is the total Mass
	Mtot = cellMass.sum()

	# loosely define the ang velocity
	# give proper size of lx
	# Its before bulk velocity subtraction
	lx = y*vz - vy*z
	ly = z*vx - x*vz
	lz = x*vy - y*vx
############################################################################
	if (ShellSphere_Bulk):
		# Find the Bulk Velocity of each shell, prior to removing from each cell in the shell
		print "Calculating the Bulk Velocity by sphere inside shell"
		ts = time.time()
		st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
		print st

		mbin = numpy.zeros(bins)
		menc = numpy.zeros(bins)
		vx_sphere_bulk_bin = numpy.zeros(bins)
		vy_sphere_bulk_bin = numpy.zeros(bins)
		vz_sphere_bulk_bin = numpy.zeros(bins)
		vx_bulkvelocity_bin = numpy.zeros(bins)
		vy_bulkvelocity_bin = numpy.zeros(bins)
		vz_bulkvelocity_bin = numpy.zeros(bins)
		
		lgradiusMin = math.log10( radiusMin)
		lgradiusSph = math.log10( radiusSphere)
		for i in range(r.size) : 
			if( r[i]/parsec < radiusMin) :
				continue
			index = int((math.log10(r[i]/parsec)-lgradiusMin)*bins/(lgradiusSph - lgradiusMin))
			if(index >= 0 and index < bins) :
				# Get the mass of each shell
				mbin[index] = mbin[index] + cellMass[i]
				
				# Find the Bulk velocity (Momentum) of each shell
				vx_sphere_bulk_bin[index] = vx_sphere_bulk_bin[index] + (vx[i]*cellMass[i])
				vy_sphere_bulk_bin[index] = vy_sphere_bulk_bin[index] + (vy[i]*cellMass[i])
				vz_sphere_bulk_bin[index] = vz_sphere_bulk_bin[index] + (vz[i]*cellMass[i])
		menc[0] = mbin[0]
		
		for shell in range(1,bins):
			menc[shell] = mbin[shell] + menc[shell-1]
			vx_sphere_bulk_bin[shell] = vx_sphere_bulk_bin[shell] + vx_sphere_bulk_bin[shell-1] 
			vy_sphere_bulk_bin[shell] = vy_sphere_bulk_bin[shell] + vy_sphere_bulk_bin[shell-1] 
			vz_sphere_bulk_bin[shell] = vz_sphere_bulk_bin[shell] + vz_sphere_bulk_bin[shell-1] 

		# Set the bulk velocity to be this bulk velocity		
		vx_bulkvelocity_bin = vx_sphere_bulk_bin/menc
		vy_bulkvelocity_bin = vy_sphere_bulk_bin/menc
		vz_bulkvelocity_bin = vz_sphere_bulk_bin/menc

###########################################################################
	if (Shell_Bulk) or (NoParticle_Bulk):
		# Find the Bulk Velocity of each shell, prior to removing from each cell in the shell
		print 'Finding Bulk Velocity by Shell'
		ts = time.time()
		st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
		print st

		mbin = numpy.zeros(bins)
		vx_shell_bulk_bin = numpy.zeros(bins)
		vy_shell_bulk_bin = numpy.zeros(bins)
		vz_shell_bulk_bin = numpy.zeros(bins)
		vx_bulkvelocity_bin = numpy.zeros(bins)
		vy_bulkvelocity_bin = numpy.zeros(bins)
		vz_bulkvelocity_bin = numpy.zeros(bins)
		
		lgradiusMin = math.log10( radiusMin)
		lgradiusSph = math.log10( radiusSphere)
		for i in range(r.size) :
			if( r[i]/parsec < radiusMin) :
				continue
			index = int((math.log10(r[i]/parsec)-lgradiusMin)*bins/(lgradiusSph - lgradiusMin))
			if(index >= 0 and index < bins) :
				# Get the mass of each shell
				mbin[index] = mbin[index] + cellMass[i]
				# Find the Bulk Velocity of each shell
				vx_shell_bulk_bin[index] = vx_shell_bulk_bin[index] + (vx[i]*cellMass[i])
				vy_shell_bulk_bin[index] = vy_shell_bulk_bin[index] + (vy[i]*cellMass[i])
				vz_shell_bulk_bin[index] = vz_shell_bulk_bin[index] + (vz[i]*cellMass[i])


       		vx_shell_bulk_bin = vx_shell_bulk_bin/mbin
	       	vy_shell_bulk_bin = vy_shell_bulk_bin/mbin
		vz_shell_bulk_bin = vz_shell_bulk_bin/mbin

		# Set the bulk velocity to be this bulk velocity
		vx_bulkvelocity_bin = vx_shell_bulk_bin
		vy_bulkvelocity_bin = vy_shell_bulk_bin
		vz_bulkvelocity_bin = vz_shell_bulk_bin


############################################################################
	print 'Obtaining vr'
	ts = time.time()
	st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
	print st

	rhobin = numpy.zeros(bins)
	volbin = numpy.zeros(bins)
	vrbin = numpy.zeros(bins)
	mdotbin = numpy.zeros(bins)

	vxbin = numpy.zeros(bins)
	vybin = numpy.zeros(bins)
	vzbin = numpy.zeros(bins)
	vr_nomass = numpy.zeros(vx.size)
	# get the radial velocity
	lgradiusMin = math.log10( radiusMin)
	lgradiusSph = math.log10( radiusSphere)
	for i in range(r.size) : 
		if( r[i]/parsec < radiusMin) :
			continue
		index = int((math.log10(r[i]/parsec)-lgradiusMin)*bins/(lgradiusSph - lgradiusMin))
		if(index >= 0 and index < bins) :
			# This is a volume weighted density. i.e. Calculate the mass
			# We'll then divide the mass by volbin
			rhobin[index] = rhobin[index] + cellMass[i]
			volbin[index] = volbin[index] + cellVolume[i]
			# reset vx, vy, vz to remove bulk velocity
			vx[i] = vx[i] - vx_bulkvelocity_bin[index]
			vy[i] = vy[i] - vy_bulkvelocity_bin[index]
			vz[i] = vz[i] - vz_bulkvelocity_bin[index]
			vr_nomass[i] = (vx[i]*x[i] + vy[i]*y[i] + vz[i]*z[i])/r[i]
			vr = vr_nomass[i] * cellMass[i]
			mdotbin[index] = mdotbin[index] + vr/r[i]
			vrbin[index] = vrbin[index] + vr
		
	vrbin = vrbin/mbin
	rhobin = rhobin/volbin

#####################################################################
	print 'Obtaining the angular momentum'
	ts = time.time()
	st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
	print st

	mbin = numpy.zeros(bins) # mbin calculated above
	menc = numpy.zeros(bins)
	angXbin = numpy.zeros(bins)
	angYbin = numpy.zeros(bins)
	angZbin = numpy.zeros(bins)
	vphi_magbin = numpy.zeros(bins)

	angX_shellsphere_bin = numpy.zeros(bins)
	angY_shellsphere_bin = numpy.zeros(bins)
	angZ_shellsphere_bin = numpy.zeros(bins)

	lx_shellsphere = numpy.zeros(bins)
	ly_shellsphere = numpy.zeros(bins)
	lz_shellsphere = numpy.zeros(bins)
	ang_velocity = numpy.zeros(bins)

	for i in range(r.size) : 
		if( r[i]/parsec < radiusMin) :
			continue
		index = int((math.log10(r[i]/parsec)-lgradiusMin)*bins/(lgradiusSph - lgradiusMin))
		if(index >= 0 and index < bins) :
			mbin[index] = mbin[index] + cellMass[i]
			# Now calculate the angular momentum (technically just r x v here)
			lx[i] = y[i]*vz[i] - vy[i]*z[i]
			ly[i] = z[i]*vx[i] - x[i]*vz[i]
			lz[i] = x[i]*vy[i] - y[i]*vx[i]
	
			angXbin[index] = lx[i] * cellMass[i] + angXbin[index]
			angYbin[index] = ly[i] * cellMass[i] + angYbin[index]
			angZbin[index] = lz[i] * cellMass[i] + angZbin[index]


	angX_sphere_bin = angXbin
	angY_sphere_bin = angYbin
	angZ_sphere_bin = angZbin

	menc[0] = mbin[0]
	for shell in range(1,bins):
		menc[shell] = mbin[shell] + menc[shell-1]
		angX_sphere_bin[shell] = angX_sphere_bin[shell] + angX_sphere_bin[shell - 1]
		angY_sphere_bin[shell] = angY_sphere_bin[shell] + angY_sphere_bin[shell - 1]
		angZ_sphere_bin[shell] = angZ_sphere_bin[shell] + angZ_sphere_bin[shell - 1]

	lx_shellsphere = angX_sphere_bin/menc
	ly_shellsphere = angY_sphere_bin/menc
	lz_shellsphere = angZ_sphere_bin/menc

	norm = numpy.sqrt(lx_shellsphere**2 + ly_shellsphere**2 + lz_shellsphere**2)

	lgrbin = lgradiusMin + (lgradiusSph-lgradiusMin)*numpy.arange(bins)/bins
	rbin = 1e1**lgrbin
	rbinpar = rbin *parsec

	# This is lmagbin (r*v)
	vphi_magbin = numpy.sqrt(lx_shellsphere**2 + ly_shellsphere**2 + lz_shellsphere**2)
	# This is v = l/r
	vphi_magbin = vphi_magbin / rbinpar

#####################################################################
	print 'Obtaining the rms velocity'
	ts = time.time()
	st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
	print st


	# get the rms velocity
	mbin = numpy.zeros(bins) # Reset mbin as its used in calculating vkep
	vrmsbin = numpy.zeros(bins)
	vrmsnbin = numpy.zeros(bins)
	vmagbin = numpy.zeros(bins)
	vmagnbin = numpy.zeros(bins)
	nbin = numpy.zeros(bins)

	for i in range(r.size) : 
		if( r[i]/parsec < radiusMin) :
			continue
		index = int((math.log10(r[i]/parsec)-lgradiusMin)*bins/(lgradiusSph - lgradiusMin))
		if(index >= 0 and index < bins) :
			mbin[index] = mbin[index] + cellMass[i]
			nbin[index] = nbin[index] + 1
			rad = r[i]
			if( rad < 1e-5) : 
				 rad = 1e-5
			# This is the correct turbulent Velocity
			vrmsn = math.sqrt((vx[i]-vrbin[index]*x[i]/rad)**2. + (vy[i]-vrbin[index]*y[i]/rad)**2. + (vz[i]-vrbin[index]*z[i]/rad)**2.)
			vrms = vrmsn*cellMass[i]
			vmagn = math.sqrt(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i])
			vmag = vmagn*cellMass[i]

			vr = vr_nomass[i]
			vrmsbin[index] = vrmsbin[index] + vrms
			vrmsnbin[index] = vrmsnbin[index] + vrmsn
			vmagbin[index] = vmagbin[index] + vmag
			vmagnbin[index] = vmagnbin[index] + vmagn

	vrmsbin = vrmsbin/mbin
	vrmsnbin = vrmsnbin/nbin
	vmagbin = vmagbin/mbin
	vmagnbin = vmagnbin/nbin

###############################################################################

	# get the Kepler velocity
	print 'Obtaining the Kepler Velocity'
	ts = time.time()
	st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
	print st

	mTbin = menc
	# mTbin is in grams
	mTbin = mTbin + particleMass * Msun # Particle Mass was in solar masses
	vKbin = numpy.sqrt( G*mTbin/rbinpar)

	numpy.savetxt(fileout, zip(rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, rhobin, mdotbin, norm, angXbin, angYbin, angZbin, vphi_magbin), fmt="%15.9E")

##############################################################################################

def ortho_grid_mapping(pf, xc, yc, zc, vxc, vyc, vzc, fileout="rad_junk.out", radiusMin=1e-3, radiusSphere=3.0, particleMass = 0.0):
	#ortho_ray(axis, coords, fields=None, pf=None, **kwargs)
	# axis: 0 = x, 1 = y, 2 = z
	# What is the Beam width yt chooses? With AMR how do I specify that then?
	# Coords appear to be in cm from 0 to 16 parsec
	plot_out_prefix = 'histogram'
	oray = pf.h.ortho_ray(2, (xc, yc))
	vz_column = oray["z-velocity"]
	dens_column = oray["Density"]
	dz = oray["dz"]
	weight_n = dens_column * dz

	mask_ = numpy.any(weight_n < 0.05)
	print mask_
	print weight_n[weight_n<0.05]
#	print weight_n
#	print dz
	#numpy.histogram(a, bins=10, range=None, normed=False, weights=None, density=None)
#	hist, npbin_edges = numpy.histogram(vz_column[weight_n < 0.05]/1e5, bins=100, weights=weight_n[weight_n < 0.05], density=True)
	hist, npbin_edges = numpy.histogram(vz_column[weight_n > 0.0005]/1e5, bins=100, weights=weight_n[weight_n > 0.0005], density=True)
	#hist, npbin_edges = numpy.histogram(np.log10(weight_n), bins=100, density=True)
	plt.bar(npbin_edges[:-1], hist, width = 0.1)
	plt.xlim(min(npbin_edges), max(npbin_edges))
#	print hist
#	print npbin_edges[:-1]

	#plt.plot(hist[0])
	plot_field = 'num'
	fileout="{0}_{1}_{2:04d}_{3}_{4:03d}.pdf".format(plot_out_prefix, quad, i, plot_field, j)
	plt.savefig(fileout)
	sys.exit()
	plt.clf()
	#plt.hist([1, 2, 1], bins=[0, 1, 2, 3])
	#hist(x, bins=10, range=None, normed=False, weights=None, cumulative=False, bottom=None, histtype=u'bar', align=u'mid', orientation=u'vertical', rwidth=None, log=False, color=None, label=None, stacked=False, hold=None, **kwargs)
	hist, npbin_edges = numpy.histogram(vz_column/1e5, bins=100, weights=None, density=True)
	plt.bar(npbin_edges[:-1], hist, width = 0.1)
	plt.xlim(min(npbin_edges), max(npbin_edges))
	#plt.plot(hist[0])
	plot_field = 'unweight'
	fileout="{0}_{1}_{2:04d}_{3}_{4:03d}.pdf".format(plot_out_prefix, quad, i, plot_field, j)
	plt.savefig(fileout)

	plt.clf()
	(mu, sigma) = scinorm.fit(vz_column/1e5)
	# the histogram of the data
	#n, bins, patches = plt.hist(datos, 60, normed=1, facecolor='green', alpha=0.75)
	n, bins, weighting = plt.hist(vz_column/1e5, bins = 100, normed =True, weights=dens_column, color='g')

	# add a 'best fit' line
	y = mlab.normpdf( bins, mu, sigma)
	#plt.plot(bins, y, 'r--', linewidth=2)
	plot_field = 'matplot'
	fileout="{0}_{1}_{2:04d}_{3}_{4:03d}.pdf".format(plot_out_prefix, quad, i, plot_field, j)
	plt.savefig(fileout)
	plt.clf()
	#numpy.savetxt(fileout, zip(vz_column, dens_column), fmt="%15.9E")
#	for i in range(r.size) : 
#		if( r[i]/parsec < radiusMin) :
#			continue

def Proj_plotting(pf, xc, yc, zc, zoom_width, quad, i, j):
	plot_out_prefix = 'movieframe'
	plot_axis = 'x'
	plot_field = plot_axis + '-velocity'
	fileout="{0}_{1}_{2:04d}_{3}_{4}_{5:03d}_{6}pc.png".format(plot_out_prefix, quad, i, plot_axis, plot_field, j, zoom_width)
        p_plot = ProjectionPlot(pf, plot_axis, plot_field, weight_field = 'Density', center = (xc, yc, zc), width = (zoom_width, 'pc'))
	#p_plot.save(fileout)
	pid = os.getpid()
        p_plot.save("junk{0:06d}".format(pid))
        shutil.move("junk{0:06d}_Projection_{1}_{2}_Density.png".format(pid, plot_axis, plot_field),fileout)

################################################

import argparse
parser = argparse.ArgumentParser(description = "start number to end number")

parser.add_argument('start', metavar='N1', type=int)
parser.add_argument('end', metavar='N2', type=int)
parser.add_argument('step', metavar='N2', type=int)
parser.add_argument('--smallsphere', action='store_true')
parser.add_argument('--bigsphere', action='store_true')
parser.add_argument('--particle', action='store_true')
parser.add_argument('--noparticle', action='store_true')
parser.add_argument('--shell', action='store_true')
parser.add_argument('--shellsphere', action='store_true')

args = parser.parse_args()

YT_CALL = False
PY_CALL = False
# These settings eventually go to yt calls
withSmallSphere = args.smallsphere
withBigSphere = args.bigsphere
withParticle = args.particle
withNoParticle = args.noparticle

Sphere_Bulk = False
Particle_Bulk = False
NoParticle_Bulk = False

if (withSmallSphere):
	compare_file = 'smallsphere'
	Sphere_Bulk = True
	Bulk_sphere_radius = 0.01
	YT_CALL = True

if (withBigSphere):
	compare_file = 'bigsphere'
	Sphere_Bulk = True
	Bulk_sphere_radius = 3.0
	YT_CALL = True

if (withParticle):
	compare_file = 'part'
	Particle_Bulk = True
	YT_CALL = True

if (withNoParticle):
	compare_file = 'nopart'
	NoParticle_Bulk = True
	YT_CALL = True

# These calls eventually go to the python written outputs
withShell = args.shell
withShellSphere = args.shellsphere

Shell_Bulk = False
ShellSphere_Bulk = False

if (withShell):
	compare_file = 'shell'
	Shell_Bulk = True
	PY_CALL = True


if (withShellSphere):
	compare_file = 'shellsphere'
	ShellSphere_Bulk = True
	PY_CALL = True

# These are universal
#plt_prefix = "BB_hdf5_chk_"
plt_prefix = "BB_hdf5_plt_cnt_"
out_prefix = "rad_profile"
quad = os.getcwd()[-5:]
# The file looks like this:
#'{out_prefix}{framestep}_{compare_file}_{particle_number}.out'                                                                                               
# i.e. rad_profile_0218_part_000.out   
# compare file options area: bigsphere, smallsphere, part, nopart, shell, shellsphere
for i in range(args.start,args.end,args.step) :	
	print 'On File: {0}_{1:04d}'.format(plt_prefix, i)
	if not (withNoParticle):
		pf = load("{0}{1:04d}".format( plt_prefix, i))

		dd = pf.h.all_data()
		xp = dd["particle_posx"]
		yp = dd["particle_posy"]
		zp = dd["particle_posz"]
		vxp = dd["particle_velocity_x"]
		vyp = dd["particle_velocity_y"]
		vzp = dd["particle_velocity_z"]
		partMass = dd["ParticleMassMsun"]

		for j in range(xp.size) :
			xc = xp[j]	 
			yc = yp[j]	 
			zc = zp[j]
			
			vxc = vxp[j]
			vyc = vyp[j]
			vzc = vzp[j]	 
			if (withParticle) or (withSmallSphere) or (withBigSphere):
				print 'Using yt'
				print 'On particle:', j + 1, 'of:', xp.size
				getRadialProfile_yt(pf,xc,yc,zc, vxc, vyc, vzc, fileout="{0}_{1}_{2:04d}_{3}_{4:03d}.out".format( out_prefix, quad, i, compare_file, j), radiusSphere=3.0, particleMass = partMass[j])

			if (withShell) or (withShellSphere):
				print 'going to python script'
				print 'On particle:', j + 1, 'of:', xp.size
				ortho_grid_mapping(pf, xc, yc, zc, vxc, vyc, vzc, fileout="test_map.out", radiusMin=1e-3, radiusSphere=3.0, particleMass = 0.0)
				#getRadialProfile_py(pf, xc, yc, zc, fileout="{0}_{1}_{2:04d}_{3}_{4:03d}.out".format( out_prefix, quad, i, compare_file, j), radiusMin=1e-3, radiusSphere=3.0, particleMass = partMass[j])
				zoom_width = 0.1
				Proj_plotting(pf, xc, yc, zc, zoom_width, quad, i, j)

	if (withNoParticle):
		# Use this for before star particle formation.
		# Make sure getRadialProfile_yt has bulk set to sphere and not particle Velocity
		pf = load("{0}{1:04d}".format( plt_prefix, i))
		dd = pf.h.all_data()
		max= dd.quantities["MaxLocation"]("Density")
		maxDens = max[0]/3e-22
		maxLoc = numpy.array(max[2:5])/3e18
		xc = max[2]
		yc = max[3]
		zc = max[4]
		vxc = 0
		vyc = 0
		vzc = 0
		#getRadialProfile_yt(pf,xc,yc,zc, vxc, vyc, vzc, fileout="{0}{1:04d}_nopart.out".format( out_prefix, i), radiusSphere=3.0)
		#getRadialProfile_yt(pf,xc,yc,zc, vxc, vyc, vzc, fileout="{0}_{1}_{2:04d}_{3}.out".format( out_prefix, quad, i, compare_file), radiusSphere=3.0, particleMass = 0.0)
		getRadialProfile_py(pf, xc, yc, zc, fileout="{0}_{1}_{2:04d}_{3}.out".format( out_prefix, quad, i, compare_file), radiusMin=1e-3, radiusSphere=3.0, particleMass = 0.0)
