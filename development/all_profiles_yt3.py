"""
This script looks at simulation data from RAMSES (FLASH suppport in a future update) and finds all star particles that have been formed.
It then creates a sphere of a given radius (default 3pc) around each star particle or density peak, and loads the information of all cells in that sphere.
It then subtracts the bulk motion velocity, options include: having yt use the partice velocity --particle, letting yt do the removal --bigsphere (radius 3pc),
our current standard --shellsphere (subtract bulk velocity by spherical lograthmic shells, --smallsphere using a sphere of radius 1e-3pc centered on particle,
After bulk velocity subraction we obtain the infall velocity vr. Then the angular momentum of the sphere.
Finally, this script then finds the vrms by taking the velocity of each cell and subtracting out the infall and angular velocity.
All data is then saved into a numpy zipped txt file, typically rad_profile*.out
"""

### All import statements ########
import matplotlib
matplotlib.use("Agg")

import yt
import matplotlib.pyplot as p

import math
import numpy
import random
import time
import datetime
import os
import numpy.linalg
import sys
import glob
from yt import YTArray
from shutil import copyfile

def getRadialProfile_yt(pf, xc, yc, zc, vxc, vyc, vzc, fileout="rad_profile.out", radiusMin=1e-3, radiusSphere=3.0, ParticleMass = 0.0) : 
	if (Bulk_by_Sphere_in_Shell): #Bigsphere
		sp = pf.h.sphere((xc, yc, zc), radiusSphere/pf['pc'])
		sp2 = pf.h.sphere((xc, yc, zc), Bulk_sphere_radius/pf['pc'])
		print Bulk_sphere_radius
		bulk = sp2.quantities["BulkVelocity"]()
	if (Bulk_by_Particle): 
		sp = pf.h.sphere((xc, yc, zc), radiusSphere/pf['pc'])
		#sp2 = pf.h.sphere((xc, yc, zc), 0.01/pf['pc'])
		bulk = numpy.array([vxc, vyc, vzc])
#	if (NoParticle_Bulk):
#		sp = pf.h.sphere((xc, yc, zc), 4.0/pf['pc'])
#		bulk = sp.quantities["BulkVelocity"]()
	sp.set_field_parameter("bulk_velocity", bulk)
	rp = BinnedProfile1D( sp, bins, "Radiuspc", 1e1**logRhoMin, radiusSphere, log_space=True)
	rp.add_fields("Density", weight="CellVolume")
	maxDens = rp["Density"][0]
	
	rp.add_fields("RadialVelocity")  
	rp.add_fields("VelocityMagnitude")
	
	rp.add_fields("CellMassMsun", accumulation=True,weight=None)

	rbin = rp["Radiuspc"]
	vrbin = rp["RadialVelocity"]
	mTbin = rp["CellMassMsun"] + ParticleMass
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

def getRadialProfile_py(pf, Particle_attributes, ParticleID, creation_time, current_time, fileout="rad_profile.out", radiusMin=1e-3, radiusSphere=3.0, ParticleMass = 0.0) : 
	global Penrose_matrix_Particles
	print 'Sphere radius is: ', radiusSphere
	ts = time.time()
	st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
	print st

	#Pull out the Particle_attributes
	xc = Particle_attributes[0]
	yc = Particle_attributes[1]
	zc = Particle_attributes[2]
	vxc = Particle_attributes[3]
	vyc = Particle_attributes[4]
	vzc = Particle_attributes[5]
	lxc = Particle_attributes[6]
	lyc = Particle_attributes[7]
	lzc = Particle_attributes[8]

	print xc, yc, zc
	print vxc, vyc, vzc
	print lxc, lyc, lzc

	sp = pf.h.sphere([xc, yc, zc], (radiusSphere, 'pc'))

	#Determine the x,y,z distances from the sink particle for each cell in the sphere.
	x = sp["x"].in_cgs() - xc
	y = sp["y"].in_cgs() - yc
	z = sp["z"].in_cgs() - zc
	r = numpy.sqrt(x*x + y*y + z*z)
	#Convert away from the YTArray
	x = numpy.array(x)
	y = numpy.array(y)
	z = numpy.array(z)
	r = numpy.array(r)

	#Calculate the l_hat of the sink particle
	lrc = numpy.sqrt(lxc*lxc + lyc*lyc + lzc*lzc)
	l_hat = (lxc, lyc, lzc)/lrc
	print 'l_hat', l_hat
#	sys.exit()
	
	if( args.magnetic) :
		Bx = sp["ramses_magx"]#In gauss
		By = sp["ramses_magy"]
		Bz = sp["ramses_magz"]
		Btot = numpy.sqrt(Bx*Bx + By*By + Bz*Bz)
		#Convert away from the YTArray
		Bx = numpy.array(Bx)
		By = numpy.array(By)
		Bz = numpy.array(Bz)
		Btot = numpy.array(Btot)
		#print 'B Field', Btot
#		sys.exit()

	# grams
	cellMass = numpy.array(sp["cell_mass"])# yt now loads in g not Msolar
#	cellMass = sp["cell_mass"]
#	cellMass = numpy.array(cellMass)
	# cgs
	dens = sp["Density"].in_cgs()#yt loads in code density
	dens = numpy.array(dens)
	# cm**3
	cellVolume = sp["cell_volume"].in_cgs()
	cellVolume = numpy.array(cellVolume)
	# These are in cm/s
	#vx = sp["x-velocity"]
	vx = sp["velocity_x"]
	vy = sp["velocity_y"]
	vz = sp["velocity_z"]
	vx = numpy.array(vx)
	vy = numpy.array(vy)
	vz = numpy.array(vz)
	
	# This is the total Mass
	Mtot = cellMass.sum()

	print 'Loosely defining ang velocity'
	ts = time.time()
	st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
	print st

	# loosely define the ang velocity
	# To give proper size of lx
	# Its before bulk velocity subtraction
	lx = y*vz - vy*z
	ly = z*vx - x*vz
	lz = x*vy - y*vx

############################################################################
	if (Bulk_by_Sphere_in_Shell): #(ShellSphere_Bulk) or (NoParticle_Bulk) or (withSmallSphere):
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
		v_vector_bulk = numpy.sqrt(vx_bulkvelocity_bin[-1]**2 + vy_bulkvelocity_bin[-1]**2 + vz_bulkvelocity_bin[-1]**2)
		print 'bulk velocity is (cm/s): ', v_vector_bulk
		print 'radius of sphere is (pc) ', radiusSphere
		print 'crossing time in sec: ', radiusSphere * parsec / v_vector_bulk


###########################################################################
###########################################################################
	if (Bulk_by_Shell):#(Shell_Bulk):
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
############################################################################
	if( args.magnetic) :
		print 'Obtaining Magnetic field etc...'
		ts = time.time()
		st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
		print st
		Bxbin = numpy.zeros(bins)
		Bybin = numpy.zeros(bins)
		Bzbin = numpy.zeros(bins)
		Btotbin = numpy.zeros(bins)

		lgradiusMin = math.log10( radiusMin)
		lgradiusSph = math.log10( radiusSphere)
		for i in range(r.size) : 
			if( r[i]/parsec < radiusMin) :
				continue
			index = int((math.log10(r[i]/parsec)-lgradiusMin)*bins/(lgradiusSph - lgradiusMin))
			if(index >= 0 and index < bins) :
				#print Bxbin[index]
				#print Bx[i]
				# To Do MASS WEIGHTED! OR NOT?
				Bxbin[index] = Bxbin[index] + Bx[i]
				Bybin[index] = Bybin[index] + By[i]
				Bzbin[index] = Bzbin[index] + Bz[i]
				Btotbin[index] = Btotbin[index] + Btot[i]
#	sys.exit()
############################################################################
############################################################################
	print 'Obtaining vr'
	ts = time.time()
	st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
	print st

	#mbin = numpy.zeros(bins) # mbin calculated above
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
			if( args.jet) :
				# x[i],y[i],z[i] are already determined as the distance from sink to that cell i
				# and r[i] is the mag of that vector distance.
				#The jet is pushing out ~1/3 of the mass.
				# Calculate the r_hat from sink particle to each cell
				r_hat = (x[i], y[i], z[i])/r[i] # cells x 3 matrix.
				dot_r_hat_l_hat = abs(r_hat[0]*l_hat[0]+r_hat[1]*l_hat[1]+r_hat[2]*l_hat[2])
				#neg_dot_r_hat_l_hat = abs(r_hat[0]*-1.0*l_hat[0]+r_hat[1]*-1.0*l_hat[1]+r_hat[2]*-1.0*l_hat[2])
				if( dot_r_hat_l_hat > jet_open_angle):
					# Inside the Jet opening angle, don't use u_r
					continue

			#mbin[index] = mbin[index] + cellMass[i]
			# This is a volume weighted density. i.e. Calculate the mass
			# We'll then divide the mass by volbin
			rhobin[index] = rhobin[index] + cellMass[i] #This is actually just a mass, not a density.
			volbin[index] = volbin[index] + cellVolume[i]
			# getting a mass weighted speed out of a velocity (Mdot * distance). 
			#vr = ((vx[i] - vx_bulkvelocity_bin[index])*x[i] + (vy[i] - vy_bulkvelocity_bin[index])*y[i] + (vz[i] - vz_bulkvelocity_bin[index])*z[i])*cellMass[i]/r[i]
			# reset vx, vy, vz to remove bulk velocity
			vx[i] = vx[i] - vx_bulkvelocity_bin[index]
			vy[i] = vy[i] - vy_bulkvelocity_bin[index]
			vz[i] = vz[i] - vz_bulkvelocity_bin[index]
				
			vr_nomass[i] = (vx[i]*x[i] + vy[i]*y[i] + vz[i]*z[i])/r[i]
			vr = vr_nomass[i] * cellMass[i]
			#mdotbin[index] = mdotbin[index] + vr # vr is mdot right now
			mdotbin[index] = mdotbin[index] + vr/r[i]
			vrbin[index] = vrbin[index] + vr

			#Original Kludge
#			if( vr < 0.):
#				mdotbin[index] = mdotbin[index] + vr/r[i]
#				vrbin[index] = vrbin[index] + vr
#
			# Original
#			mdotbin[index] = mdotbin[index] + vr/r[i]
#			vrbin[index] = vrbin[index] + vr
	vrbin = vrbin/mbin
	# Check to see if these come out the same:
	#rhobin = mbin/volbin
	rhobin = rhobin/volbin#This converts it properly into a density.
	
	# Find the middle radius of the bin
	lgrbin = lgradiusMin + (lgradiusSph-lgradiusMin)*numpy.arange(bins)/bins
	rbin = 1e1**lgrbin
	# put r in cm
	rbinparsec = rbin * parsec

#####################################################################
	print 'Obtaining the angular momentum'
	ts = time.time()
	st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
	print st

	mbin = numpy.zeros(bins)
	angXbin = numpy.zeros(bins)
	angYbin = numpy.zeros(bins)
	angZbin = numpy.zeros(bins)
#	Number_of_cells_bin = numpy.zeros(bins)
#	old_vphi_magbin = numpy.zeros(bins)

	for i in range(r.size) : 
		if( r[i]/parsec < radiusMin) :
			continue
		index = int((math.log10(r[i]/parsec)-lgradiusMin)*bins/(lgradiusSph - lgradiusMin))
		if(index >= 0 and index < bins) :
			if( args.jet) :
				# x[i],y[i],z[i] are already determined as the distance from sink to that cell i
				# and r[i] is the mag of that vector distance.
				#The jet is pushing out ~1/3 of the mass.
				# Calculate the r_hat from sink particle to each cell
				r_hat = (x[i], y[i], z[i])/r[i] # cells x 3 matrix.
				dot_r_hat_l_hat = abs(r_hat[0]*l_hat[0]+r_hat[1]*l_hat[1]+r_hat[2]*l_hat[2])
				#neg_dot_r_hat_l_hat = abs(r_hat[0]*-1.0*l_hat[0]+r_hat[1]*-1.0*l_hat[1]+r_hat[2]*-1.0*l_hat[2])
				if( dot_r_hat_l_hat > jet_open_angle):
					# Inside the Jet opening angle, don't use u_r
					continue
			mbin[index] = mbin[index] + cellMass[i]
			# Now calculate the angular momentum of cell i (technically just r x v here)
			lx[i] = y[i]*vz[i] - vy[i]*z[i]
			ly[i] = z[i]*vx[i] - x[i]*vz[i]
			lz[i] = x[i]*vy[i] - y[i]*vx[i]
			#Number_of_cells_bin[index] = Number_of_cells_bin[index] + 1
		# Don't need this as we have reset vx vy and vz to include this
			#lx[i] = y[i]*(vz[i] - vz_bulkvelocity_bin[index]) - (vy[i] - vy_bulkvelocity_bin[index])*z[i]
			#ly[i] = z[i]*(vx[i] - vx_bulkvelocity_bin[index]) - (vz[i] - vz_bulkvelocity_bin[index])*x[i]
			#lz[i] = x[i]*(vy[i] - vy_bulkvelocity_bin[index]) - (vx[i] - vx_bulkvelocity_bin[index])*y[i]
			# vphi is perp to this ang momentum vector
			angXbin[index] = lx[i] * cellMass[i] + angXbin[index]
			angYbin[index] = ly[i] * cellMass[i] + angYbin[index]
			angZbin[index] = lz[i] * cellMass[i] + angZbin[index]

	# Set the Angular momentum vectors before getting specific ang momentum
	Lbin = numpy.matrix([angXbin, angYbin, angZbin])
	# This is the mass weighted specific ang momentum per bin not by sphere
	angXbin = angXbin/mbin
	angYbin = angYbin/mbin
	angZbin = angZbin/mbin
#	L_mag = angZbin
#	L_mag = numpy.sqrt(angXbin**2 + angYbin**2 + angZbin**2)
	norm = numpy.sqrt(angXbin**2 + angYbin**2 + angZbin**2)
#	old_vphi_magbin = numpy.sqrt(angXbin**2 + angYbin**2 + angZbin**2)
#	old_vphi_magbin = old_vphi_magbin / rbinparsec
########################################################
	print 'Calculating Moment of Inertia'
	ts = time.time()
	st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
	print st

	Ixxbin = numpy.zeros(bins)
	Iyybin = numpy.zeros(bins)
	Izzbin = numpy.zeros(bins)
	Ixybin = numpy.zeros(bins)
	Ixzbin = numpy.zeros(bins)
	Iyxbin = numpy.zeros(bins)
	Iyzbin = numpy.zeros(bins)
	Izxbin = numpy.zeros(bins)
	Izybin = numpy.zeros(bins)
	#I_matrixbin = numpy.zeros(bins)
	Omega = numpy.empty((bins), dtype=numpy.object)

	# If using speedup.so comment out this for loop as this is done in fortran.
	for i in range(r.size) : 
		if( r[i]/parsec < radiusMin) :
			continue
		index = int((math.log10(r[i]/parsec)-lgradiusMin)*bins/(lgradiusSph - lgradiusMin))
		if(index >= 0 and index < bins) :
			if( args.jet) :
				# x[i],y[i],z[i] are already determined as the distance from sink to that cell i
				# and r[i] is the mag of that vector distance.
				#The jet is pushing out ~1/3 of the mass.
				# Calculate the r_hat from sink particle to each cell
				r_hat = (x[i], y[i], z[i])/r[i] # cells x 3 matrix.
				dot_r_hat_l_hat = abs(r_hat[0]*l_hat[0]+r_hat[1]*l_hat[1]+r_hat[2]*l_hat[2])
				#neg_dot_r_hat_l_hat = abs(r_hat[0]*-1.0*l_hat[0]+r_hat[1]*-1.0*l_hat[1]+r_hat[2]*-1.0*l_hat[2])
				if( dot_r_hat_l_hat > jet_open_angle):
					# Inside the Jet opening angle, don't use u_r
					continue

			Ixx = (y[i]*y[i] + z[i]*z[i]) * cellMass[i]
			Iyy = (x[i]*x[i] + z[i]*z[i]) * cellMass[i]
			Izz = (x[i]*x[i] + y[i]*y[i]) * cellMass[i]
			Ixy = -x[i]*y[i] * cellMass[i]
			Ixz = -x[i]*z[i] * cellMass[i]
				#Iyx = Ixy
			Iyz = -y[i]*z[i] * cellMass[i]
				#Izx = Ixz
				#Izy = Iyz
			Ixxbin[index] = Ixxbin[index] + Ixx
			Iyybin[index] = Iyybin[index] + Iyy
			Izzbin[index] = Izzbin[index] + Izz
			Ixybin[index] = Ixybin[index] + Ixy
			Ixzbin[index] = Ixzbin[index] + Ixz
			Iyzbin[index] = Iyzbin[index] + Iyz
#	import speedup
#	Ixxbin, Iyybin, Izzbin, Ixybin, Ixzbin, Izybin = speedup.moment_of_inertia( r,x,y,z,cellMass,parsec,radiusMin,lgradiusMin,lgradiusSph,bins)
	# Set these outside the for loop
	Iyxbin = Ixybin
	Izxbin = Ixzbin
	Izybin = Iyzbin
	print "Finished calculating the moment of Inertia elements"
	ts = time.time()
	st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
	print st

#	test_dot_product = 0

	for index in range(bins):
		# Pick out Matrix elements for this bin
		Ixxbin_index =  Ixxbin[index]
		Ixybin_index =  Ixybin[index]
		Ixzbin_index =  Ixzbin[index]

		Iyxbin_index =  Iyxbin[index]
		Iyybin_index =  Iyybin[index]
		Iyzbin_index =  Iyzbin[index]

		Izxbin_index =  Izxbin[index]
		Izybin_index =  Izybin[index]
		Izzbin_index =  Izzbin[index]
		
		# Collapse all matrix elements of this bin into I for this bin
		I_matrixbin = numpy.matrix([[Ixxbin_index, Ixybin_index, Ixzbin_index], [Iyxbin_index, Iyybin_index, Iyzbin_index], [Izxbin_index, Izybin_index, Izzbin_index]])
		# This line is necessary for the no particle  plt files
		# In those cases, we have a 3x3 matrix with 0 values,
		# and Invert I has a hard time with that.
		if Ixxbin_index == Ixybin_index == Ixzbin_index == Iyybin_index == Iyzbin_index == Izzbin_index == 0:
			Omega_X = 0
			Omega_Y = 0
			Omega_Z = 0
			Omega[index] = numpy.matrix([[Omega_X],[ Omega_Y], [Omega_Z]])
			continue
		# Invert I
		#print I_matrixbin
		# Need to place a try except block in here, as
		# occasionally we see a singular matrix.
		try:
			I_invert_matrix = numpy.linalg.inv(I_matrixbin)
		except numpy.linalg.LinAlgError:
			print "Linear algebra error, can't invert this matrix."
			print "Will try using Penrose-Moore psuedo-inverse to find a 'best' solution,"
			print "But you should probably avoid this particle"
#			Penrose_matrix_Particles = Penrose_matrix_Particles.append(ParticleID)
			try:
				I_invert_matrix = numpy.linalg.pinv(I_matrixbin)
			except numpy.linalg.LinAlgError:#numpy.linalg.linalg.LinAlgError:
				print "Failed. Code should hard crash here."
#				pass_to_std_out("Failed Inversions of I. Code should hard crash here.")
#				I_invert_matrix = I_matrixbin*0.0
#		I_invert_matrix = numpy.linalg.inv(I_matrixbin)			
		# L
		L_this_bin = Lbin[0,index], Lbin[1,index], Lbin[2,index]
		# This gives Omega = I_inverse * L for this bin
		Omega_this_bin = numpy.dot(I_invert_matrix, L_this_bin)
		Omega_X = Omega_this_bin[0, 0]
		Omega_Y = Omega_this_bin[0, 1]
		Omega_Z = Omega_this_bin[0, 2]

		Omega[index] = numpy.matrix([[Omega_X],[ Omega_Y], [Omega_Z]])

#####################################################################
	print 'Obtaining the rms velocity'
	ts = time.time()
	st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
	print st


	# get the rms velocity
	mbin = numpy.zeros(bins) # Reset mbin as its used in calculating vkep
	vrmsbin = numpy.zeros(bins)
	vrms_r_bin = numpy.zeros(bins)
	vrms_l_bin = numpy.zeros(bins)
	vrmsnbin = numpy.zeros(bins)
	vmagbin = numpy.zeros(bins)
	vmagnbin = numpy.zeros(bins)
	nbin = numpy.zeros(bins)
	vphi_magbin = numpy.zeros(bins)
	mvphi_magbin = numpy.zeros(bins)
	vrms_phi_bin = numpy.zeros(bins)
	vrms_theta_bin = numpy.zeros(bins)


	for i in range(r.size) : 
		if( r[i]/parsec < radiusMin) :
			continue
		index = int((math.log10(r[i]/parsec)-lgradiusMin)*bins/(lgradiusSph - lgradiusMin))
		#index = int(r[i]*bins/(radiusSphere*parsec))
		if(index >= 0 and index < bins) :
			if( args.jet) :
				# x[i],y[i],z[i] are already determined as the distance from sink to that cell i
				# and r[i] is the mag of that vector distance.
				#The jet is pushing out ~1/3 of the mass.
				# Calculate the r_hat from sink particle to each cell
				r_hat = (x[i], y[i], z[i])/r[i] # cells x 3 matrix.
				dot_r_hat_l_hat = abs(r_hat[0]*l_hat[0]+r_hat[1]*l_hat[1]+r_hat[2]*l_hat[2])
				#neg_dot_r_hat_l_hat = abs(r_hat[0]*-1.0*l_hat[0]+r_hat[1]*-1.0*l_hat[1]+r_hat[2]*-1.0*l_hat[2])
				if( dot_r_hat_l_hat > jet_open_angle):
					# Inside the Jet opening angle, don't use u_r
					continue

			mbin[index] = mbin[index] + cellMass[i]
			nbin[index] = nbin[index] + 1
			rad = r[i]
			if( rad < 1e-5) : 
				rad = 1e-5
			# Pull out the terms for Omega for this shell
			Omega_X = Omega[index][0]
			Omega_Y = Omega[index][1]
			Omega_Z = Omega[index][2]
			# Convert them to floats from numpy.array objects
			Omega_X = float(Omega_X)
			Omega_Y = float(Omega_Y)
			Omega_Z = float(Omega_Z)
			# Calculate V_omega for x y and z
			v_Omega_X = (Omega_Y * z[i] - Omega_Z * y[i])
			v_Omega_Y = -(Omega_X * z[i] - Omega_Z * x[i])
			v_Omega_Z = (Omega_X * y[i] - Omega_Y * x[i])
			# Calculate the magnitude of V_omega to compare to vphi_magbin
			v_Omega = numpy.sqrt(v_Omega_X**2 + v_Omega_Y**2 + v_Omega_Z**2)
			#print v_Omega / 1e5
			mvphi_magbin[index] = mvphi_magbin[index] + v_Omega * cellMass[i]

			# This is the correct turbulent Velocity
			#vrmsn = math.sqrt((vx[i]-vrbin[index]*x[i]/rad)**2. + (vy[i]-vrbin[index]*y[i]/rad)**2. + (vz[i]-vrbin[index]*z[i]/rad)**2.)
			vrmsn = math.sqrt((vx[i]-v_Omega_X-vrbin[index]*x[i]/rad)**2. + (vy[i]-v_Omega_Y-vrbin[index]*y[i]/rad)**2. + (vz[i]-v_Omega_Z-vrbin[index]*z[i]/rad)**2.)
			vrms_x = vx[i]-v_Omega_X-vrbin[index]*x[i]/rad
			vrms_y = vy[i]-v_Omega_Y-vrbin[index]*y[i]/rad
			vrms_z = vz[i]-v_Omega_Z-vrbin[index]*z[i]/rad
			vrms = vrmsn*cellMass[i]
			vmagn = math.sqrt(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i])
			vmag = vmagn*cellMass[i]
			# Calculate the radial portion of the turbulent velocity.
			vrms_radial_x = vrms_x * x[i]/rad 
			vrms_radial_y = vrms_y * y[i]/rad
			vrms_radial_z = vrms_z * z[i]/rad
			#vrms_radial = numpy.sqrt(vrms_radial_x**2 + vrms_radial_y**2 + vrms_radial_z**2) * cellMass[i]
			vrms_radial = (vrms_radial_x + vrms_radial_y + vrms_radial_z)

			#sin_theta_ = numpy.sqrt(rad**2 - z**2)/rad
			vrms_theta_x = vrms_x * (x[i]*z[i])/(rad*numpy.sqrt(rad**2 - z[i]**2))
			vrms_theta_y = vrms_y * (y[i]*z[i])/(rad*numpy.sqrt(rad**2 - z[i]**2))
			vrms_theta_z = vrms_z * -1.0 * numpy.sqrt(rad**2 - z[i]**2)/rad
			#vrms_theta = numpy.sqrt(vrms_theta_x**2 + vrms_theta_y**2 + vrms_theta_z**2) * cellMass[i]
			vrms_theta = (vrms_theta_x + vrms_theta_y + vrms_theta_z)

			vrms_phi_x = vrms_x * -1.0 * y[i]/numpy.sqrt(rad**2 - z[i]**2)
			vrms_phi_y = vrms_y * x[i]/(numpy.sqrt(rad**2 - z[i]**2))
			vrms_phi_z = 0.0
			#vrms_phi = numpy.sqrt(vrms_phi_x**2 + vrms_phi_y**2 + vrms_phi_z**2) * cellMass[i]
			# Can't sum in quad?
			vrms_phi = (vrms_phi_x + vrms_phi_y + vrms_phi_z)
			# avg out using the abs for phi and theta
			#vrms_phi = (vrms_phi)
			#vrms_lateral = numpy.sqrt(vrms_theta**2 + vrms_phi**2)

			vrmsbin[index] = vrmsbin[index] + vrms
			vrmsnbin[index] = vrmsnbin[index] + vrmsn
			vmagbin[index] = vmagbin[index] + vmag
			vmagnbin[index] = vmagnbin[index] + vmagn

			vrms_r_bin[index] = vrms_r_bin[index] + vrms_radial**2  * cellMass[i]
			#vrms_l_bin[index] = vrms_l_bin[index] + vrms_lateral
			vrms_theta_bin[index] = vrms_theta_bin[index] + vrms_theta**2  * cellMass[i]
			vrms_phi_bin[index] = vrms_phi_bin[index] + vrms_phi**2  * cellMass[i]

#	print 'vrms', vrmsbin
#	print 'mbin', mbin
#	sys.exit()
	vrmsbin = vrmsbin/mbin
	vrmsnbin = vrmsnbin/nbin
	vmagbin = vmagbin/mbin
	vmagnbin = vmagnbin/nbin
	vphi_magbin = mvphi_magbin / mbin

	vrms_r_bin = numpy.sqrt(vrms_r_bin/mbin)
	#vrms_l_bin = vrms_l_bin/mbin
	vrms_theta_bin = numpy.sqrt(vrms_theta_bin/mbin)
	vrms_phi_bin = numpy.sqrt(vrms_phi_bin/mbin)
	vrms_l_bin = (vrms_theta_bin + vrms_phi_bin) / 2.0
###############################################################################

	# get the Kepler velocity
	print 'Obtaining the Kepler Velocity'
	ts = time.time()
	st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
	print st

	menc = numpy.zeros(bins)
	menc[0] = mbin[0]
	for shell in range(1,bins):
		menc[shell] = mbin[shell] + menc[shell-1]
	mTbin = menc
	# mTbin is in grams
	mTbin = mTbin + ParticleMass * Msun # Particle Mass was in solar masses
	# calculated rbin in the ang momentum section
	# and converted from being in parsecs
	#lgrbin = lgradiusMin + (lgradiusSph-lgradiusMin)*numpy.arange(bins)/bins
	#rbin = 1e1**lgrbin
	#vKbin = numpy.sqrt( G*mTbin*Msun/((rbin)*parsec))
	#vKbin = numpy.sqrt( G*mTbin/((rbin)*parsec))
	vKbin = numpy.sqrt( G*mTbin/rbinparsec)

	sum_velocities = numpy.sqrt(vphi_magbin**2 + vrbin**2 + vrmsbin**2)
	
	# Include the particle mass
	Particle_mass = numpy.full((bins), ParticleMass)
	CreationTime = numpy.full((bins), creation_time)
	current_time = numpy.full((bins), current_time)
	if( args.magnetic):
		print fileout[:-4] + 'magnetic.out'
		numpy.savetxt(fileout[:-4] + 'magnetic.out', zip(Bxbin, Bybin, Bzbin, Btotbin), fmt="%15.9E")
	numpy.savetxt(fileout, zip(rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, rhobin, mdotbin, norm, angXbin, angYbin, angZbin, vphi_magbin, sum_velocities, Particle_mass, CreationTime, current_time, vrms_r_bin, vrms_l_bin, vrms_theta_bin, vrms_phi_bin), fmt="%15.9E")

################################################
################################################
################################################
################################################
################################################

def getDensityExclude( pf, xc, yc, zc, xexclude, yexclude, zexclude, fileoutrho_PDF="rho_PDF.out", radiusExclude=1.0, exclude=True) : 
	sp = pf.h.all_data()

	dens = sp["Density"]
	cellMass = sp["CellMassMsun"]
	x = sp["x"]
	y = sp["y"]
	z = sp["z"]
	
	# create a pdf
	volWeightedPDF = numpy.zeros(bins)
	massWeightedPDF = numpy.zeros(bins)
	volWeightedIPDF = numpy.zeros(bins)
	massWeightedIPDF = numpy.zeros(bins)

	rE2 = (radiusExclude*parsec)**2.

	for i in range(dens.size) :
		includeDens = True
		excludeDens = False
		if( exclude ) :
			for j in range(xexclude.size) :
				xpos = x[i] - xexclude[j]
				ypos = y[i] - yexclude[j]
				zpos = z[i] - zexclude[j]
				r2 = xpos*xpos + ypos*ypos + zpos*zpos
				if( r2 < rE2) :
					includeDens = False
					excludeDens = True
					break
		if( includeDens) : 	
			rho = dens[i]
			mcell = cellMass[i]
			lRho = math.log10(rho/meanDens)
			index = int((lRho - logRhoMin)*bins/(logRhoMax-logRhoMin))
			if(index >= 0 and index < bins) :
				volWeightedPDF[index] = volWeightedPDF[index] + 1
				massWeightedPDF[index] = massWeightedPDF[index] + mcell
		if( excludeDens) : 	
			rho = dens[i]
			mcell = cellMass[i]
			lRho = math.log10(rho/meanDens)
			index = int((lRho - logRhoMin)*bins/(logRhoMax-logRhoMin))
			if(index >= 0 and index < bins) :
				volWeightedIPDF[index] = volWeightedIPDF[index] + 1
				massWeightedIPDF[index] = massWeightedIPDF[index] + mcell
	# normalized 
	volWeightedPDF = volWeightedPDF/volWeightedPDF.sum()
	massWeightedPDF = massWeightedPDF/massWeightedPDF.sum()
	volWeightedIPDF = volWeightedIPDF/volWeightedIPDF.sum()
	massWeightedIPDF = massWeightedIPDF/massWeightedIPDF.sum()

	# write out
	lRhoArr = numpy.arange(bins)*(logRhoMax - logRhoMin)/bins + logRhoMin
	numpy.savetxt(fileoutrho_PDF, zip( lRhoArr, volWeightedPDF, massWeightedPDF,volWeightedIPDF, massWeightedIPDF), fmt="%15.9E")

def getDensityFunctions( pf, xc, yc, zc, fileoutrho_r="rho_r.out", fileoutrho_PDF="rho_PDF.out", radiusSphere=3.0, allData = False) : 
	sp = 0
	if( allData ) : 
		sp = pf.h.all_data()
	else :
		sp = pf.h.sphere([xc,yc,zc], radiusSphere/pf['pc'])

	dens = sp["Density"]
	cellMass = sp["CellMassMsun"]
	
	# create a pdf
	volWeightedPDF = numpy.zeros(bins)
	massWeightedPDF = numpy.zeros(bins)
	
	for i in range(dens.size) :
		rho = dens[i]
		mcell = cellMass[i]
		lRho = math.log10(rho/meanDens)
		index = int((lRho - logRhoMin)*bins/(logRhoMax-logRhoMin))
		if(index >= 0 and index < bins) :
			volWeightedPDF[index] = volWeightedPDF[index] + 1
			massWeightedPDF[index] = massWeightedPDF[index] + mcell
	
	# normalized 
	volWeightedPDF = volWeightedPDF/volWeightedPDF.sum()
	massWeightedPDF = massWeightedPDF/massWeightedPDF.sum()

	# write out
	lRhoArr = numpy.arange(bins)*(logRhoMax - logRhoMin)/bins + logRhoMin
	numpy.savetxt(fileoutrho_PDF, zip( lRhoArr, volWeightedPDF, massWeightedPDF), fmt="%15.9E")

	if( allData) :
		return

	# calculate rho as function of r
	x = sp["x"] - xc
	y = sp["y"] - yc
	z = sp["z"] - zc
	r = numpy.sqrt(x*x + y*y + z*z)

	volWeightedRho_r = numpy.zeros(bins)
	massWeightedRho_r = numpy.zeros(bins)
	vol_r = numpy.zeros(bins)
	mass_r = numpy.zeros(bins)
	
	for i in range(dens.size) :
		rho = dens[i]
		mcell = cellMass[i]
		dist = r[i] 
		index = int(bins*dist/(radiusSphere*parsec))
		if( index >= 0 and index < bins) : 
			volWeightedRho_r[index] = volWeightedRho_r[index] + rho
			massWeightedRho_r[index] = massWeightedRho_r[index] + rho*mcell
			vol_r[index] = vol_r[index] + 1
			mass_r[index] = mass_r[index] + mcell
	
	volWeightedRho_r = volWeightedRho_r/vol_r
	massWeightedRho_r = massWeightedRho_r/mass_r

	# write out
	rbin = radiusSphere*numpy.arange(bins)/bins
	numpy.savetxt(fileoutrho_r, zip( rbin, volWeightedRho_r, massWeightedRho_r), fmt="%15.9E")

def getLarsonsLaw(pf, xc, yc, zc, fileout="larson.out", radiusSphere=3.0,trials=10,randPoints = True) : 
	sp = pf.h.sphere([xc,yc,zc], radiusSphere/pf['pc'])

	x = sp["x"] - xc
	y = sp["y"] - yc
	z = sp["z"] - zc

	vx = sp["x-velocity"]
	vy = sp["y-velocity"]
	vz = sp["z-velocity"]
	vrmsbin = numpy.zeros(bins)
	nbin = numpy.zeros(bins)

	if ( not randPoints) : 
		trials = 1

	for trial in range(trials) : 
		ridx = 0
		if(not randPoints) :
			for i in range(x.size) :
				r2 = (x[i]*x[i] + y[i]*y[i] + z[i]*z[i])/(parsec*parsec)
				if ( r2 < 0.03*0.03) :
					ridx = i 
				
		else : 
			ridx = random.randint(0,x.size-1) 

		vxr = vx[ridx]
		vyr = vy[ridx]
		vzr = vz[ridx]
		xr = x[ridx]
		yr = y[ridx]
		zr = z[ridx]

		for i in range(x.size) : 
			if(i != ridx) : 
				dist = math.sqrt((x[i] - xr)**2. + (y[i] - yr)**2 + (z[i] - zr)**2)
				vrms = (vx[i] - vxr)**2. + (vy[i] - vyr)**2 + (vz[i] - vzr)**2

				index = int(bins*dist/(2.*radiusSphere*parsec))

				if(index < bins) :
					vrmsbin[index] = vrmsbin[index] + vrms
					nbin[index] = nbin[index] + 1
		

	vrmsbin = numpy.sqrt(vrmsbin/nbin)
	rbin = 2.*radiusSphere*numpy.arange(bins)/bins

	numpy.savetxt(fileout, zip(rbin, vrmsbin), fmt="%12.7e")

########################################################################################

def FLASH_file_existance_chk(file_value):
	fn_plt = plt_prefix+"{0:04d}".format(file_value)
	fn_part = part_prefix+"{0:04d}".format(file_value)
	file_plt_exist = glob.glob(fn_plt)
	file_part_exist = glob.glob(fn_part)
	if not file_plt_exist:
		print 'File: "', fn_plt, '" does not exist, moving to the next file'
		return False
	elif not file_part_exist:
		print 'File: "', fn_part, '" does not exist, moving to the next file'
		return False
	else:
		return True

def FLASH_additional_stuff():
	#fileout="{0}_{1}_{2:04d}_{3}_{4}.out".format( out_prefix, quad, int(File_number), compare_file, int(ParticleID))
	if (withMaxDensity):
		# Use this for before star particle formation.
				#pf = yt.load("{0}{1:04d}".format( plt_prefix, i))
		dd = pf.h.all_data()
		max= dd.quantities["MaxLocation"]("Density")
		maxDens = max[0]/3e-22
		maxLoc = numpy.array(max[2:5])/3e18
		xc = max[2]
		yc = max[3]
		zc = max[4]
		vxc = 0
#			vyc = 0
#			vzc = 0
#			current_time = pf.current_time
#			ParticleID = int(42)
#			creation_time = 0.0
#			Particle_age = 0.0
#			this_Particle_Mass = 0.0
#			print 'Passing to the analysis script'
#			pass_to_analysis_script(pf, xc, yc, zc, ParticleID, creation_time, 
#						current_time, radiusMin, radiusSphere, ParticleMass, 
#						File_number, sinkfile, file)


def FLASH_load_particle_info(plt_file):
	#FLASH
	pf = yt.load(plt_file)
	dd = pf.h.all_data()
	xp = dd["particle_posx"]
	yp = dd["particle_posy"]
	zp = dd["particle_posz"]
	vxp = dd["particle_velocity_x"]
	vyp = dd["particle_velocity_y"]
	vzp = dd["particle_velocity_z"]
	partMass = dd["ParticleMassMsun"]
	creation_time = dd["particle_creation_time"]
	current_time = pf.current_time
	particle_ID_list = dd["particle_index"]

def FLASH_obtain_individual_part_attributes(part_lst_num, File_number, pf, Particle_ID_list, xp, yp, zp, vxp, vyp, vzp, partMass, creation_time, current_time):
	xc = xp[part_lst_num]	 
	yc = yp[part_lst_num]	 
	zc = zp[part_lst_num]
	vxc = vxp[part_lst_num]
	vyc = vyp[part_lst_num]
	vzc = vzp[part_lst_num]
	ParticleID = Particle_ID_list[part_lst_num]
	ParticleMass = partMass[part_lst_num]
	#print ParticleID
	this_creation_time = creation_time[part_lst_num]
	Particle_age = current_time - this_creation_time
	Particle_age = Particle_age / (numpy.pi*1e7)
	return xc, yc, zc, vxc, vyc, vzc, ParticleID, this_creation_time, ParticleMass

#############################################
#############################################
#############################################
#############################################
#############################################


def print_to_text_file(txt_filename, File_number, ParticleID, xc, yc, zc, creation_time, current_time, ParticleMass):
	# Fill with the new values for this Particle and timestep.
	File_number = numpy.full((1), File_number)
	ParticleID = numpy.full((1), ParticleID)
	ParticleMass = numpy.full((1), ParticleMass)
	xc = numpy.full((1), xc)
	yc = numpy.full((1), yc)
	zc = numpy.full((1), zc)
	creation_time = numpy.full((1), creation_time)
	current_time = numpy.full((1), current_time)
	with open(txt_filename, "a") as Particle_location:
		numpy.savetxt(Particle_location, zip(File_number, ParticleID, ParticleMass, xc, yc, zc, creation_time, current_time))

def standardize_File_number(input_File_number):
	output_File_number = "%05d"%input_File_number
	if len(output_File_number) >= 6:
		print 'This system assumes less than one million output files.'
		sys.exit()
	return output_File_number

def RAMSES_obtain_individual_part_attributes(pf, index, Particle_ID_list, partMass, xstar, ystar, zstar, vxstar, vystar, vzstar, lxstar, lystar, lzstar):
	ParticleID = Particle_ID_list[index]
	ParticleMass = partMass[index]
	xc = pf.quan(xstar[index], "cm")
	yc = pf.quan(ystar[index], "cm")
	zc = pf.quan(zstar[index], "cm")
	vxc = vxstar[index]
	vyc = vystar[index]
	vzc = vzstar[index]
	lxc = lxstar[index]
	lyc = lystar[index]
	lzc = lzstar[index]
	return ParticleID, ParticleMass, xc, yc, zc, vxc, vyc, vzc, lxc, lyc, lzc

def Obtain_Particles(num_var_packed, pf, sinkfile, File_number, sinkcsvfile):#, xc, yc, zc):#, Init_Restart):
	global Init_Restart
	global xc_search
	global yc_search
	global zc_search
	#Init for this particle & timestep
	Particle_ID_list = numpy.float64(0.0)
	partMass = 0.0
	r_star = 0.0
	poly_n = 0.0
	md = 0.0
	polystate = 0
	pjet = 0.0
	xstar = 0.0
	ystar = 0.0
	zstar = 0.0
	vxstar = 0.0
	vystar = 0.0
	vzstar = 0.0
	lxstar = 0.0
	lystar = 0.0
	lzstar = 0.0
	Particle_age = 0.0
#	print 'obtaining particles', len( num_var_packed)
	if len( num_var_packed) == 0 or ( args.point):
#	if True:	
		dd = pf.all_data()
		if ( args.maxdensity):
			# If there is no particles, but we chose to look for the highest density point
			# This sets the required variables that we would obtain from the sinkfile
			# It then continues as if we have just a single particle in this file.
			print 'Finding max density location.'
			max= dd.quantities["MaxLocation"]("Density")
		elif ( args.backwards):
			#This is the case if we are tracking a particle backwards
			Particle_ID_list = numpy.float64(withParticleIDValue)
			if ( withRestart and Init_Restart):
				#Restarted the run, read from particle location file.
				#ld_File_number, ld_PartID, ld_PartMass, ld_xc, ld_yc, ld_zc, ld_creation, ld_current
				ld_File_number, ld_xc, ld_yc, ld_zc = numpy.loadtxt("{0}/particle_{1}_location.txt".format(output_location, args.ParticleID), usecols=[0,3,4,5], unpack=True)
				for value in range(len(ld_File_number)):
					ld_file = standardize_File_number(ld_File_number[value])
					if ( ld_file == File_number):
						print ld_file, File_number
						xc_search = ld_xc[value]
						yc_search = ld_yc[value]
						zc_search = ld_zc[value]
						Init_Restart = False
						break
			print 'Backwards search location: xc, yc, zc:', xc_search, yc_search, zc_search
			sp = pf.sphere(YTArray( [xc_search, yc_search, zc_search], "cm"), (0.5, 'pc'))
			max = sp.quantities["MaxLocation"]("Density")
		else:
			Particle_ID_list = numpy.float64(withParticleIDValue)
			print 'looking at hardcoded specified point for testing.'
			xc_search = 2.38813134766e+19
			yc_search = 2.91307177734e+19
			zc_search = 1.19757666016e+19
			sp = pf.sphere(YTArray( [xc_search, yc_search, zc_search], "cm"), (0.5, 'pc'))
			max = sp.quantities["MaxLocation"]("Density")
		maxDens = pf.quan(max[0], "g*cm**(-3)")
		max_Location = numpy.array(max[1:4]) * 16. * parsec# convert to cm 3e18
		xstar = max_Location[0]
		ystar = max_Location[1]
		zstar = max_Location[2]
		print 'These are the Max Density coordinates: ', xstar, ystar, zstar
	elif len( num_var_packed) == 8: # Magnetic Field run
		Particle_ID_list, partMass, xstar, ystar, zstar, \
		    vxstar, vystar, vzstar = numpy.loadtxt( sinkfile, unpack=True, skiprows=3, comments="=")
	elif len( num_var_packed) == 9: # The Nakano run. # or no jets
		Particle_ID_list, partMass, r_star, xstar, ystar, zstar, \
		    vxstar, vystar, vzstar = numpy.loadtxt( sinkfile, unpack=True, skiprows=3, comments="=")
	elif len( num_var_packed) == 12:# The Offner run, no tracking pjet. 
		Particle_ID_list, partMass, r_star, xstar, ystar, zstar, \
		    vxstar, vystar, vzstar, poly_n, md, polystate = numpy.loadtxt( sinkfile, unpack=True, skiprows=3, comments="=")
	elif len( num_var_packed) == 13:# Offner with jet momentum.
		Particle_ID_list, partMass, r_star, xstar, ystar, zstar, \
		    vxstar, vystar, vzstar, poly_n, md, polystate, pjet = numpy.loadtxt( sinkfile, unpack=True, skiprows=3, comments="=")
	else :
		print "It appears that this file does not correspond to a Nakano run, or either version of the Offner runs."
		print len(num_var_packed)
		print num_var_packed
		sys.exit()
	if len( num_var_packed) != 0:
		# Obtain particle angular momentum from the sink.csv file.
		# Eventually, we will write this info out into sink_.info as well.
		# Note that lx, ly and lzstar are in CODE UNITS! particle_age has been converted to s.
		lxstar, lystar, lzstar, Particle_age = Obtain_Part_info_csv(sinkcsvfile)

	return Particle_ID_list, partMass, r_star, xstar, ystar, zstar, vxstar, vystar, vzstar, poly_n, md, polystate, pjet, lxstar, lystar, lzstar, Particle_age

def Obtain_Part_info_csv(sinkcsvfile):
# From Ramses output_sink.f90 for teh csv file
#     write(ilun,'(I10,12(A1,ES20.10))')idsink(isink),',',msink(isink),&
#          ',',xsink(isink,1),',',xsink(isink,2),',',xsink(isink,3),&
#          ',',vsink(isink,1),',',vsink(isink,2),',',vsink(isink,3),&
#          ',',lsink(isink,1),',',lsink(isink,2),',',lsink(isink,3),&
#          ',',t-tsink(isink),',',dMBHoverdt(isink)
	# Note that all of these values are returned in CODE UNITS! not cgs/solar like the sinkfile
	scale_t      =  0.387201000000000E+04
	Particle_ID_list, partMass, xstar, ystar, zstar, \
	    vxstar, vystar, vzstar, lxstar, lystar, lzstar, \
	    Particle_Age, dMBhoverdt_star = numpy.loadtxt( sinkcsvfile, unpack=True, skiprows=0, delimiter=',', comments="=")
	Particle_Age = Particle_Age * scale_t
	return	lxstar, lystar, lzstar, Particle_Age

###########################################################
########_____START_OF_Particle_Reduction_____########
###########################################################

def Particle_Reduction(index):
	# First, check that the required files exist
	#FLASH
	if ( withFLASH4):
		if not (FLASH_file_existance_chk):
			return
		print 'On File: {0}{1:04d}'.format(plt_prefix, index)
		# Flash Stuff Here
		plt_file = '{0}{1:04d}'.format(plt_prefix, index)
		FLASH_load_Particle_info(plt_file)

	else:
		#RAMSES
		print 'On Folder: ' + cwd + '/output_{0:05d}'.format(index)
		# Create the 5 digit file number so that we can reference it for plotting
		File_number = standardize_File_number(index)
		infofile = prefix + "{0}/info_{0}.txt".format(File_number)
		sinkfile = prefix + "{0}/sink_{0}.info".format(File_number)
		sinkcsvfile = prefix + "{0}/sink_{0}.csv".format(File_number)
		if( not os.path.isfile( sinkfile)) or ( not os.path.isfile( infofile)) or ( not os.path.isfile( sinkcsvfile)):
			print "Either the sink_.info, info_.txt or sink_.csv file is missing."
			return

		# Copy infofile and sinkfile to the output location, they are used in the plotting routine.
		copyfile(infofile, '{0}/info_{1}.txt'.format(output_location, File_number))
		copyfile(sinkfile, '{0}/sink_{1}.info'.format(output_location, File_number))
		copyfile(sinkcsvfile, '{0}/sink_{1}.csv'.format(output_location, File_number))
		#The files exist, time to load up the information
		pf = yt.load(infofile)
		num_var_packed = numpy.loadtxt( sinkfile, unpack=True, skiprows=3, comments="=")
		if len( num_var_packed) == 0:
			# If there is no particles in this timestep
			#and we have not chosen to look for the highest density point jump out.
			# not withMaxdensity is true, so it doesn't matter about the trackback command.
			if not ( args.backwards):
				if not ( withMaxDensity):
					print 'No particle, jumping out.'
					return
			else:
				Particle_existance = False
		# Obtain particle information from the sinkfile
		# The current total possible options to be packed in these ramses pieces thus far.
		if( args.magnetic) :
			print("Length unit: ", pf.length_unit)
			print("Time unit: ", pf.time_unit)
			print("Mass unit: ", pf.mass_unit)
			print("Velocity unit: ", pf.velocity_unit)
			pf.add_field(("gas", "ramses_magx"), function=_ramses_magx, units="gauss") 
			pf.add_field(("gas", "ramses_magy"), function=_ramses_magy, units="gauss") 
			pf.add_field(("gas", "ramses_magz"), function=_ramses_magz, units="gauss") 
			print pf.all_data()["ramses_magx"]

		Particle_ID_list, partMass, r_star, xstar, ystar, zstar, \
		    vxstar, vystar, vzstar, poly_n, md, polystate, pjet, \
		    lxstar, lystar, lzstar, Particle_age = Obtain_Particles(num_var_packed, pf, sinkfile, File_number, sinkcsvfile)
#		print lxstar, lystar, lzstar
#		print Particle_age
#		sys.exit()
		# Now that we've loaded the relevant data, pull out what we require.
		#Creation time is available in FLASH, not implemented yet for the RAMSES data.
		creation_time = 0.0
		current_time = pf.current_time
		#Particle_age = 0.0 # current_time - this_creation_time

		# Now that we have all the information associated with all particles in this timestep,
		# How are we requested to loop through?
		if( args.maxdensity):
			Particle_list = 1
		elif xstar.size == 1:# If there are zero particles, or one, both trigger here.
			Particle_list = 1
		elif xstar.size >= 2:
			half_point = xstar.size/2
			if ( args.first) or ( args.second):
				print 'will only do half of the list of particles.'
				Particle_list = xstar.size/2 + 1
			elif ( withParallel):
				#TO DO Figure out how this works and make it work.
				print 'Running the script in parallel.'
				parallel_size = int(math.ceil(1.0 * xstar.size/size)) 
				Particle_list = parallel_size - rank
				#for j in range(rank, parallel_size, size):
				#	print 'On particle:', j + 1, 'of:', xstar.size, ' in File: {0}{1:04d}'.format(plt_prefix, index)
				#	if j  >= xstar.size:
				#		break
			else: # Do all on this cpu
				Particle_list = xstar.size
		for j in range(Particle_list):
			if Particle_list == 1:
				if Particle_ID_list == args.ParticleID or (withAllParticles):
					# This converts to the proper data types needed.
					ParticleID = int(Particle_ID_list)
					ParticleMass = partMass
					xc = pf.quan(xstar, "cm")
					yc = pf.quan(ystar, "cm")
					zc = pf.quan(zstar, "cm")
					vxc = vxstar
					vyc = vystar
					vzc = vzstar
					lxc = lxstar
					lyc = lystar
					lzc = lzstar
				else:
					print 'The ID of the only particle in this timestep does not match the requested ID.'
					print 'PartID', str(int(Particle_ID_list)) + ', ', 'Requested ID', args.ParticleID
					continue
			else:
				#For all and args.first, we start at j = 0
				if ( args.second):
					j = j+half_point
				elif ( withParallel):
					#TO DO Figure out the proper stepping.
					j = j + size
				if j + 1 > xstar.size:
					# This would query past the end length of xstar.size.
					break
				if Particle_ID_list[j] == withParticleIDValue or (withAllParticles):
					ParticleID, ParticleMass, xc, yc, zc, vxc, vyc, vzc, lxc, lyc, lzc = \
					    RAMSES_obtain_individual_part_attributes(pf, j, Particle_ID_list, \
					    partMass, xstar, ystar, zstar, vxstar, vystar, vzstar, lxstar, lystar, lzstar)
				else:
					continue

			if ( args.backwards):
				if ( ParticleID == args.ParticleID):
					# Save the information for the current timestep and particle.
					# If we need to restart when prior to the particles formation, this will allow us
					# to jump to the location of its density peak.
					txt_filename = '{0}/particle_{1}_location.txt'.format(output_location, int(ParticleID))
					print_to_text_file(txt_filename, File_number, ParticleID, xc, yc, zc, creation_time, current_time, ParticleMass)
					# If we are tracing backwards update the search coordinates
					# To center on the particle's current location
					global xc_search
					global yc_search
					global zc_search
					xc_search = xc
					yc_search = yc
					zc_search = zc
			print 'Particle Coordinates', xc, yc, zc
			#Package up the Particle Attributes to pass to Radial Profiles.
			Particle_attributes = [xc, yc, zc, vxc, vyc, vzc, lxc, lyc, lzc]
			fileout="{0}/{1}_{2}_{3}_{4}.out".format(output_location, out_prefix, File_number, compare_file, int(ParticleID))
			# Finally, pass the appropriate data to the analysis script.
			if ("particle" in compare_file):# or ("bigsphere" in compare_file):
				getRadialProfile_yt(pf,xc,yc,zc, vxc, vyc, vzc, fileout, radiusSphere, ParticleMass)
			else:
				getRadialProfile_py(pf, Particle_attributes, ParticleID, creation_time, current_time, fileout, radiusMin, radiusSphere, ParticleMass)
			if ( withFLASH4):
				print 'Finished Analysing particle:', j + 1, 'of:', xstar.size, 'in File:', cwd, '/{0}{1:04d}'.format(plt_prefix, index)
			else:
				print 'Finished Analysing particle:', j + 1, 'of:', xstar.size, 'in Folder:', cwd, '/output_{0:05d}'.format(index)
	# After looping through all particles in this time step.
	return


def _ramses_magx(field, data):
	return yt.YTArray((data["x-Bfield-left"] + data["x-Bfield-right"]) * numpy.sqrt(4.*numpy.pi) * 0.5 / 3872.01, 'gauss')

def _ramses_magy(field, data):
	return yt.YTArray((data["y-Bfield-left"] + data["y-Bfield-right"]) * numpy.sqrt(4.*numpy.pi) * 0.5 / 3872.01, 'gauss')

def _ramses_magz(field, data):
	return yt.YTArray((data["z-Bfield-left"] + data["z-Bfield-right"]) * numpy.sqrt(4.*numpy.pi) * 0.5 / 3872.01, 'gauss')

#def _ramses_magy(field, data, pf):
#	return yt.YTArray((data["y-Bfield-left"] + data["y-Bfield-right"]) * numpy.sqrt(4.*numpy.pi) * 0.5 / 3872.01, 'gauss')
#
#def _ramses_magz(field, data, pf):
#	return yt.YTArray((data["z-Bfield-left"] + data["z-Bfield-right"]) * numpy.sqrt(4.*numpy.pi) * 0.5 / 3872.01, 'gauss')
#



###########################################################
###########################################################
###########################################################
import argparse
parser = argparse.ArgumentParser(description = "start number to end number, stride length, reduction method")
parser.add_argument('start', metavar='N1', type=int, help ='Start value for hdf5 files')
parser.add_argument('end', metavar='N2', type=int, help='End Value for hdf5 files, note runs until End-1')
parser.add_argument('step', metavar='N3', type=int, help='Stride length')
parser.add_argument('ParticleID', metavar='N4', type=int, nargs='?', default=0, 
		    help='Particle ID you want to reduce, use 0 for no particle density peak.')
parser.add_argument('bulk_vel_method', metavar='N5', type=str, nargs='?', default='shellsphere',
		    help='method of removing the bulk motion of the gas. Options are: shellsphere, bigsphere, smallsphere, particle, shell.')

parser.add_argument('--maxdensity', action='store_true')
parser.add_argument('--chk', action='store_true')
parser.add_argument('--allparticles', action='store_true')
parser.add_argument('--jet', action='store_true')
parser.add_argument('--backwards', action='store_true')
parser.add_argument('--restart', action='store_true')
parser.add_argument('--first', action='store_true')
parser.add_argument('--second', action='store_true')
parser.add_argument('--parallel', action='store_true')
parser.add_argument('--FLASH4', action='store_true')
parser.add_argument('--magnetic', action='store_true')
parser.add_argument('--point', action='store_true', help='This allows you to look at a specific location. coords are hardcoded.')

args = parser.parse_args()
withParticleIDValue = args.ParticleID
withMaxDensity = args.maxdensity
withCheckpoint = args.chk
withAllParticles = args.allparticles
withTrackBackwards = args.backwards
withRestart = args.restart
withParallel = args.parallel
withFLASH4 = args.FLASH4

#########Inits####################

bins = 70
meanDens = 3e-22
i = 0
Msun = 1.99e33
G = 6.67e-8
parsec = 3.09e18
logRhoMin = -3.0

alpha_jet = 1.0 # Fudge factor for the opening angle of the jet range it 0.9 to 1.1
# theta0_jet=0.3d0 # From RAMSES jet_parameters.f90 #Note radians
theta_jet = 0.3
jet_open_angle = math.cos( alpha_jet * theta_jet) #math.cos expects radians.
print 'jet open angle', jet_open_angle
#sys.exit()
#Depending on if we are reducing the data on scinet or another location, choose the output file path accordingly
cwd = os.getcwd()
output_location = cwd + '/python_output'

prefix = "output_"
out_prefix = "rad_profile"
# The default behaviour is setup for shellsphere.
radiusMin = 1e-3 # In parsecs
radiusSphere = 3.0 # In parsecs
Bulk_by_Sphere_in_Shell = False
Bulk_by_Shell = False
Bulk_by_Particle = False # This is used solely in yt
bulk_vel_accepted_strings = {'shellsphere', 'bigsphere', 'smallsphere', 'particle', 'shell'}

# This list is for any I matricies that are singular,
# or have issues with numpy.linalg.inv
# If they work with the penrose psuedo invert, I'd like to know
global Penrose_matrix_Particles
Penrose_matrix_Particles = []

# The output file looks like this:
#'{out_prefix}{framestep}_{compare_file}_{particle_number}.out'
# i.e. rad_profile_0218_shellsphere_000.out   

if (withParallel):
	from mpi4py import MPI
	comm = MPI.COMM_WORLD
	size = comm.Get_size()
	rank = comm.Get_rank()

#Check that the python_output folder exists, if not, create it.
if not (os.path.isdir(output_location)) :
	try:
		os.mkdir(output_location)
	except:
		print "Do not have permission to create python_output folder in current working directory"
		print cwd
		sys.exit()

if args.bulk_vel_method in bulk_vel_accepted_strings:
	compare_file = args.bulk_vel_method
	if ("shellsphere" == compare_file):
		Bulk_by_Sphere_in_Shell = True
	# These are modifications from the default shellsphere settings.
	elif ("bigsphere" == compare_file):
		radiusSphere = 5.0
	elif ('particle' == compare_file): #This goes to yt
		Bulk_by_Particle = True
	elif ("smallsphere" == compare_file):
		radiusSphere = 0.15
	elif ("shell" == compare_file):
		Bulk_by_Shell = True
else:
	print "You have not selected a valid method to remove the bulk velocity."
	print "Possible methods include: ", bulk_vel_accepted_strings
	sys.exit()


if ( withFLASH4):
	# If we wish to use the checkpoint files, rather than the output files.
	if ( withCheckpoint):
		plt_prefix = "BB_hdf5_chk_"
		part_prefix = "BB_hdf5_part_"
	else:
		plt_prefix = "BB_hdf5_plt_cnt_"
		part_prefix = "BB_hdf5_part_"

##########################################
############## Main Program ##############
##########################################
if (withTrackBackwards):
	#This is to pass xc, yc, zc back and forth with obtain particle
	#In case we wish to track a particle backwards in time, it is useful
	#to search for max density around the previous timesteps coords.
	#These are are overwritten in every case.
	global xc_search
	global yc_search
	global zc_search
	xc_search = 0.0
	yc_search = 0.0
	zc_search = 0.0
	if ( withRestart):
		global Init_Restart
		Init_Restart = True
	for i in range(args.end-1, args.start-1, -args.step) :
		print 'looping'
		Particle_Reduction(i)
		print "finished reduction on ", i
else: #Step forward in time and reduce particles
	for i in range(args.start,args.end,args.step) :	
		Particle_Reduction(i)

#print 'These particles needed a pseudo-inversion, not sure if you should trust them'
#print Penrose_matrix_particles
