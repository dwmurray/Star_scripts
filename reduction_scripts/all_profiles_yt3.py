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

#from yt.mods import * #depricated in yt3
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


bins = 70
meanDens = 3e-22
i = 0
Msun = 1.99e33
G = 6.67e-8
parsec = 3.09e18
logRhoMin = -3.0

def pass_to_std_out(output_text) :
	temp = sys.stdout #store original stdout object for later
	if withFirstHalf:
		txt_filename = '{0}/output_prints_{1}_{2}_{3}.txt'.format(filesys_location, args.start, args.end, 'first')
	elif withSecondHalf:
		txt_filename = '{0}/output_prints_{1}_{2}_{3}.txt'.format(filesys_location, args.start, args.end, 'second')
	elif withBackTracking:
		txt_filename = '{0}/output_prints_{1}_{2}_{3}.txt'.format(filesys_location, args.start, args.end, 'backwards')
	else:
		txt_filename = '{0}/output_prints_{1}_{2}.txt'.format(filesys_location, args.start, args.end)
	sys.stdout = open(txt_filename,'a') #redirect all prints to this log file
	print(output_text) #nothing appears at interactive prompt
	sys.stdout.close() #ordinary file object
	sys.stdout = temp #restore print commands to interactive prompt

def Obtain_particles(sinkfile) :
	if( not os.path.isfile( sinkfile)) :
		return
	try : 
		#ID, mass, xstar, ystar, zstar = numpy.loadtxt( sinkfile, usecols=[0,1,2,3,4], unpack=True, skiprows=3, comments="=")
		ID, mass, xstar, ystar, zstar, vxstar, vystar, vzstar = numpy.loadtxt( sinkfile, usecols=[0,1,2,3,4,5,6,7], unpack=True, skiprows=3, comments="=")
	except ValueError : 
		print 'no particles', 1+1
		pass_to_std_out('no particles')
#		ID=0
		return ID, mass, xstar, ystar, zstar, vxstar, vystar, vzstar
	#return ID, mass, xstar, ystar, zstar
	return ID, mass, xstar, ystar, zstar, vxstar, vystar, vzstar

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

def getRadialProfile_py(pf, xc, yc, zc, this_particle_id, this_creation_time, current_time, fileout="rad_profile.out", radiusMin=1e-3, radiusSphere=3.0, particleMass = 0.0) : 
	txt_2_print = 'Sphere radius is: ', str(radiusSphere)
	pass_to_std_out(txt_2_print)
	print 'Sphere radius is: ', radiusSphere
	ts = time.time()
	st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
	print st

	sp = pf.h.sphere([xc, yc, zc], (radiusSphere, 'pc'))

	x = sp["x"].in_cgs() - xc
	y = sp["y"].in_cgs() - yc
	z = sp["z"].in_cgs() - zc
	r = numpy.sqrt(x*x + y*y + z*z)
	#Convert away from the YTArray
	x = numpy.array(x)
	y = numpy.array(y)
	z = numpy.array(z)
	r = numpy.array(r)

	# grams
	cellMass = sp["cell_mass"]# * Msun #yt now loads in g not Msolar
	cellMass = numpy.array(cellMass)
	# cgs
	dens = sp["Density"].in_cgs()#yt loads in code density
	dens = numpy.array(dens)
	# cm**3
	cellVolume = sp["cell_volume"].in_cgs()# Otherwise yt loads in code density
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
	# give proper size of lx
	# Its before bulk velocity subtraction
	lx = y*vz - vy*z
	ly = z*vx - x*vz
	lz = x*vy - y*vx

############################################################################
	if (ShellSphere_Bulk) or (NoParticle_Bulk) or (withSmallSphere):
		# Find the Bulk Velocity of each shell, prior to removing from each cell in the shell
		print "Calculating the Bulk Velocity by sphere inside shell"
		ts = time.time()
		st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
		print st

		txt_2_print = "Calculating the Bulk Velocity by sphere inside shell"
		pass_to_std_out(txt_2_print)
		pass_to_std_out(st)

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
	if (Shell_Bulk):
		# Find the Bulk Velocity of each shell, prior to removing from each cell in the shell
		print 'Finding Bulk Velocity by Shell'
		ts = time.time()
		st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
		print st

		txt_2_print = "Finding Bulk Velocity by Shell"
		pass_to_std_out(txt_2_print)
		pass_to_std_out(st)

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
	txt_2_print = 'Obtaining vr'
	pass_to_std_out(txt_2_print)
	pass_to_std_out(st)
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
			#mbin[index] = mbin[index] + cellMass[i]
			# This is a volume weighted density. i.e. Calculate the mass
			# We'll then divide the mass by volbin
			rhobin[index] = rhobin[index] + cellMass[i]
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
		
	vrbin = vrbin/mbin
	# Check to see if these come out the same:
	#rhobin = mbin/volbin
	rhobin = rhobin/volbin
	
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
	txt_2_print = 'Obtaining angular momentum'
	pass_to_std_out(txt_2_print)
	pass_to_std_out(st)

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
			mbin[index] = mbin[index] + cellMass[i]
			# Now calculate the angular momentum (technically just r x v here)
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

	txt_2_print = 'Calculating Moment of Inertia.'
	pass_to_std_out(txt_2_print)
	pass_to_std_out(st)


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
			Ixx = (y[i]*y[i] + z[i]*z[i]) * cellMass[i]
			Iyy = (x[i]*x[i] + z[i]*z[i]) * cellMass[i]
			Izz = (x[i]*x[i] + y[i]*y[i]) * cellMass[i]
			Ixy = -x[i]*y[i] * cellMass[i]
			Ixz = -x[i]*z[i] * cellMass[i]
			Iyx = Ixy
			Iyz = -y[i]*z[i] * cellMass[i]
			Izx = Ixz
			Izy = Iyz
			Ixxbin[index] = Ixxbin[index] + Ixx
			Iyybin[index] = Iyybin[index] + Iyy
			Izzbin[index] = Izzbin[index] + Izz

			Ixybin[index] = Ixybin[index] + Ixy
			#Iyxbin[index] = Ixybin[index]

			Ixzbin[index] = Ixzbin[index] + Ixz
			#Izxbin[index] = Ixybin[index]

			Izybin[index] = Izybin[index] + Izy
			#Iyzbin[index] = Ixybin[index]
#	import speedup
#	Ixxbin, Iyybin, Izzbin, Ixybin, Ixzbin, Izybin = speedup.moment_of_inertia( r,x,y,z,cellMass,parsec,radiusMin,lgradiusMin,lgradiusSph,bins)
	# Set these outside the for loop
	Iyxbin = Ixybin
	Izxbin = Ixzbin
	Iyzbin = Izybin
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
			#Penrose_matrix_particles = Penrose_matrix_particles.append(this_particle_id)
			try:
				I_invert_matrix = numpy.linalg.pinv(I_matrixbin)
			except numpy.linalg.LinAlgError:#numpy.linalg.linalg.LinAlgError:
				print "Failed. Code should hard crash here."
				pass_to_std_out("Failed Inversions of I. Code should hard crash here.")
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

	txt_2_print = "Obtained Kepler Velocity"
	pass_to_std_out(txt_2_print)
	pass_to_std_out(st)

	menc = numpy.zeros(bins)
	menc[0] = mbin[0]
	for shell in range(1,bins):
		menc[shell] = mbin[shell] + menc[shell-1]
	mTbin = menc
	# mTbin is in grams
	mTbin = mTbin + particleMass * Msun # Particle Mass was in solar masses
	# calculated rbin in the ang momentum section
	# and converted from being in parsecs
	#lgrbin = lgradiusMin + (lgradiusSph-lgradiusMin)*numpy.arange(bins)/bins
	#rbin = 1e1**lgrbin
	#vKbin = numpy.sqrt( G*mTbin*Msun/((rbin)*parsec))
	#vKbin = numpy.sqrt( G*mTbin/((rbin)*parsec))
	vKbin = numpy.sqrt( G*mTbin/rbinparsec)

	sum_velocities = numpy.sqrt(vphi_magbin**2 + vrbin**2 + vrmsbin**2)
	
	# Include the particle mass
	#print particleMass
	store_particle_mass = numpy.zeros(bins)
	store_particle_mass.fill(particleMass)

	store_part_creation_time = numpy.zeros(bins)
	store_part_creation_time.fill(this_creation_time)

	store_current_time = numpy.zeros(bins)
	store_current_time.fill(current_time)

	numpy.savetxt(fileout, zip(rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, rhobin, mdotbin, norm, angXbin, angYbin, angZbin, vphi_magbin, sum_velocities, store_particle_mass, store_part_creation_time, store_current_time, vrms_r_bin, vrms_l_bin, vrms_theta_bin, vrms_phi_bin), fmt="%15.9E")

	txt_2_print = "Saved output file."
	pass_to_std_out(txt_2_print)


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


def print_to_text_file(txt_filename, four_digit_file_number, withParticleIDValue, xc, yc, zc, this_creation_time, current_time, particleMass):
#	# Reset the values to be stored to zero.
	zip_store_file_number = numpy.zeros(1)
	zip_store_particle_ID = numpy.zeros(1)
	zip_store_xc = numpy.zeros(1)
	zip_store_yc = numpy.zeros(1)
	zip_store_zc = numpy.zeros(1)
	zip_store_creation_time = numpy.zeros(1)
	zip_store_particle_mass = numpy.zeros(1)
	zip_store_current_time = numpy.zeros(1)
	# Fill with the new values for this particle and timestep.
	zip_store_file_number.fill(str(four_digit_file_number))
	zip_store_particle_ID.fill(withParticleIDValue)
	zip_store_xc.fill(xc)
	zip_store_yc.fill(yc)
	zip_store_zc.fill(zc)
	zip_store_creation_time.fill(this_creation_time)
	zip_store_particle_mass.fill(particleMass)
	zip_store_current_time.fill(current_time)
	print 'Now storing location, mass and time and particle id values.'
	#with open("particle_location.txt", "a") as particle_location:
	with open(txt_filename, "a") as particle_location:
		numpy.savetxt(particle_location, zip(zip_store_file_number, zip_store_particle_ID, zip_store_xc, zip_store_yc, zip_store_zc, zip_store_creation_time, zip_store_particle_mass, zip_store_current_time))

def standardize_file_number(input_file_number):
	if len(str(input_file_number)) == 2:
		output_file_number = str('000') + str(input_file_number)
	elif len(str(input_file_number)) == 3:
		output_file_number = str('00') + str(input_file_number)
	elif len(str(input_file_number)) == 4:
		output_file_number = str('0') + str(input_file_number)
	elif len(str(input_file_number)) == 5:
		output_file_number = str(input_file_number)
	else:
		print 'We"ve gone above 999k files apparently.'
		sys.exit()
	return output_file_number

def FLASH_load_particle_info(plt_file):
	#FLASH
	#pf = yt.load("{0}{1:04d}".format( plt_prefix, i))
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
	particle_ID = dd["particle_index"]


def flash_file_existance_chk(file_value):
	fn_plt = plt_prefix+"{0:04d}".format(file_value)
	fn_part = part_prefix+"{0:04d}".format(value)
	file_plt_exist = glob.glob(fn_plt)
	if not file_plt_exist:
		print 'File: "', fn_plt, '" does not exist, moving to the next file'
		return False #continue
	file_part_exist = glob.glob(fn_part)
	if not file_part_exist:
		print 'File: "', fn_part, '" does not exist, moving to the next file'
		return False #continue
	return True

def flash_obtain_individual_part_attributes(part_lst_num, file_number, pf, particle_ID, xp, yp, zp, vxp, vyp, vzp, partMass, creation_time, current_time):
	xc = xp[part_lst_num]	 
	yc = yp[part_lst_num]	 
	zc = zp[part_lst_num]
	vxc = vxp[part_lst_num]
	vyc = vyp[part_lst_num]
	vzc = vzp[part_lst_num]
	this_particle_id = particle_ID[part_lst_num]
	particleMass = partMass[part_lst_num]
	#print this_particle_id
	this_creation_time = creation_time[part_lst_num]
	particle_age = current_time - this_creation_time
	particle_age = particle_age / (numpy.pi*1e7)
	return xc, yc, zc, vxc, vyc, vzc, this_particle_id, this_creation_time, particleMass


def RAMSES_obtain_individual_part_attributes(part_lst_num, file_number, pf, particle_ID, xstar, ystar, zstar, vxstar, vystar, vzstar, partMass, creation_time, current_time):
	print particle_ID[part_lst_num], withParticleIDValue
	if particle_ID[part_lst_num] == withParticleIDValue or (withAllParticles):
		xc = xstar[part_lst_num]	 
		yc = ystar[part_lst_num]	 
		zc = zstar[part_lst_num]
		print xc, yc, zc, partMass[part_lst_num]
		xc = pf.quan(xc, "cm")
		yc = pf.quan(yc, "cm")
		zc = pf.quan(zc, "cm")
		vxc = vxstar[part_lst_num]
		vyc = vystar[part_lst_num]
		vzc = vzstar[part_lst_num]
		print xc, yc, zc, partMass[part_lst_num]
		particleMass = partMass[part_lst_num]
		this_particle_id = particle_ID[part_lst_num]
		print this_particle_id
		#this_creation_time = creation_time[part_lst_num]
		this_creation_time = 0.0
		particle_age = 0.0
		current_time = 0.0
		#particle_age = current_time - this_creation_time
		#particle_age = particle_age / (numpy.pi*1e7)
		print 'On particle:', part_lst_num + 1, 'of:', xstar.size, 'in File: {0}{1}'.format(prefix, file_number)
		txt_2_print = 'On particle:', part_lst_num + 1, 'of:', xstar.size, 'in File: {0}{1}'.format(prefix, file_number)
		pass_to_std_out(txt_2_print)
	else:
		print 'This particle ID does not match, ' + \
		    'may not be in this octant.'
		print 'Check the particle_location.txt file for ' + \
		    'a list of particles that are in this quadrant.'
		pass_to_std_out("PartID doesn't match. Check part location file.")
		sys.exit()
		#continue
	return xc, yc, zc, vxc, vyc, vzc, this_particle_id, this_creation_time, particleMass



def pass_to_analysis_script(pf, xc, yc, zc, this_particle_id, this_creation_time, current_time, radiusMin, radiusSphere, particleMass, file_number, sinkfile, filein):
	print 'Writing to particle location files.'
	#txt_filename = 'particle_location.txt'
	#print_to_text_file(txt_filename, file_number, this_particle_id, xc, yc, zc, this_creation_time, current_time, particleMass)
	txt_filename = '{0}/particle_{1}_location.txt'.format(filesys_location, withParticleIDValue)
	print_to_text_file(txt_filename, file_number, this_particle_id, xc, yc, zc, this_creation_time, current_time, particleMass)
	copyfile(sinkfile, '{0}/sink_{1}.info'.format(filesys_location, file_number))
	copyfile(filein, '{0}/info_{1}.txt'.format(filesys_location, file_number))

	#fileout="{0}_{1}_{2:04d}_{3}_{4}.out".format( out_prefix, quad, int(file_number), compare_file, int(this_particle_id))
	fileout="{0}/{1}_{2:04d}_{3}_{4}.out".format(filesys_location, out_prefix, int(file_number), compare_file, int(this_particle_id))
	print fileout
	if (withParticle) or (withBigSphere):
		print 'Using yt'
		getRadialProfile_yt(pf,xc,yc,zc, vxc, vyc, vzc, fileout, radiusSphere, particleMass)

	if (withShell) or (withSmallSphere) or (withShellSphere):
		print 'going to python script'
		getRadialProfile_py(pf, xc, yc, zc, this_particle_id, this_creation_time, current_time, fileout, radiusMin, radiusSphere, particleMass)


#############################################################
#########_______MAIN_TRACING_BACKWARDS_IN_TIME_______########
#############################################################
def trackbackwards():
	First_non_particle = True
	# Need end-1 and start-1 otherwise will not follow same protocol as looping forwards
	for i in range(args.end-1, args.start-1, -args.step) :
		# First, check that the two files exist:
		if withFLASH4:
			if not (flash_file_existance_chk):
				continue
			print 'On File: {0}{1:04d}'.format(plt_prefix, i)
		else:
			print 'On File: output_{0:04d}'.format(i)
		# Create the 5 digit file number so that we can reference it for plotting
		file_number = standardize_file_number(i)
		#print standardize_file_number
		if withRestart:
			print 'Restart Flag has been seen, Moving to grab particle values from particle location list.'
			particle_exists_here = False
		else:
			try:
				if withFLASH4:
					#Flash.
#					dd = pf.h.all_data()
#					xp = dd["particle_posx"]
#					yp = dd["particle_posy"]
#					zp = dd["particle_posz"]
#					vxp = dd["particle_velocity_x"]
#					vyp = dd["particle_velocity_y"]
#					vzp = dd["particle_velocity_z"]
#					partMass = dd["ParticleMassMsun"]
#					creation_time = dd["particle_creation_time"]
#					current_time = pf.current_time
#					particle_ID = dd["particle_index"]
					plt_file = '{0}{1:04d}'.format(plt_prefix, i)
					FLASH_load_particle_info(plt_file)
					pass
				else:
					#pf = load("{0}{1:04d}".format( plt_prefix, i))
					file = prefix+"{0:05d}/info_{0:05d}.txt".format(i)
					print file
					sinkfile = prefix+"{0:05d}/sink_{0:05d}.info".format(i)
					print sinkfile
					pf = yt.load(file)
					#particle_ID, partMass, xstar, ystar, zstar = Obtain_particles(sinkfile) #Ramses
					particle_ID, partMass, xstar, ystar, zstar, vxstar, vystar, vzstar = Obtain_particles(sinkfile) #Ramses
				particle_exists_here = True
#				if (particle_exists_here):
#					txt_output_file = '/home/m/murray/dwmurray/scratch/ramses_jeans_{0}/particle_location.txt'.format(filesys_location)
#					with open(txt_output_file, "a") as particle_location:
#						particle_location.write("Particle ID's in this file: " + str(particle_ID))
#						particle_location.write('\n')
#						particle_location.close()
				if xstar.size == 1:
					if particle_ID == withParticleIDValue or (withAllParticles):
						xc = xstar	 
						yc = ystar	 
						zc = zstar
						print xc, yc, zc, partMass
						xc = pf.quan(xc, "cm")
						yc = pf.quan(yc, "cm")
						zc = pf.quan(zc, "cm")
						print xc, yc, zc, partMass
						particleMass = partMass
						this_particle_id = particle_ID
						print this_particle_id
						#this_creation_time = creation_time[j]
						this_creation_time = 0.0
						particle_age = 0.0
						current_time = 0.0
					else:
						continue
					
				if xstar.size >= 2:
					for j in range(xstar.size) :
						#print particle_ID[j], withParticleIDValue
						if particle_ID[j] == withParticleIDValue or (withAllParticles):
							xc = xstar[j]	 
							yc = ystar[j]	 
							zc = zstar[j]
							#print xc, yc, zc, partMass[j]
							xc = pf.quan(xc, "cm")
							yc = pf.quan(yc, "cm")
							zc = pf.quan(zc, "cm")
						#vxc = vxstar[j]
						#vyc = vystar[j]
						#vzc = vzstar[j]
							print xc, yc, zc, partMass[j]
							particleMass = partMass[j]
							this_particle_id = particle_ID[j]
							print this_particle_id
						#this_creation_time = creation_time[j]
							h_this_part = j
							this_creation_time = 0.0
							particle_age = 0.0
							current_time = 0.0
						else:

							continue
			except KeyError:
				particle_exists_here = False
			except UnboundLocalError:
				particle_exists_here = False
		if (particle_exists_here):
#			with open("particle_location.txt", "a") as particle_location:
#				particle_location.write("Particle ID's in this file: " + str(particle_ID))
#				particle_location.write('\n')
#				particle_location.close()
#			print withAllParticles
#			for j in range(xp.size) :
#				print particle_ID[j], withParticleIDValue
#				if particle_ID[j] == withParticleIDValue or (withAllParticles):
#					xc = xp[j]	 
#					yc = yp[j]	 
#					zc = zp[j]
#					vxc = vxp[j]
#					vyc = vyp[j]
#					vzc = vzp[j]
#					this_particle_id = particle_ID[j]
#					particleMass = partMass[j]
#					this_creation_time = creation_time[j]
#					particle_age = current_time - this_creation_time
#					particle_age = particle_age / (numpy.pi*1e7)
#				else:
##					print 'This particle ID does not match, ' + \
##					    'may not be in this octant.'
##					print 'Check the particle_location.txt file for ' + \
##					    'a list of particles that are in this quadrant.'
#					#print particle_ID
#					continue
				print 'passing to analysis'
				if xstar.size >= 2:
					print 'On particle:', h_this_part + 1, 'of:', xstar.size, ' in File: {0}{1:05d}'.format(plt_prefix, i)
				print 'Particle Mass is', particleMass
				pass_to_analysis_script(pf, xc, yc, zc, this_particle_id, this_creation_time, current_time, radiusMin, radiusSphere, particleMass, file_number, sinkfile, file)

#############################################################
#########_______Continuing_BACKWARDS_IN_TIME_______########
#############################################################
		First_non_particle = True
			# This closes the particle loop inside particles exists
		if not (particle_exists_here):
			if withRestart:
				loaded_file_number, loaded_ParticleIDValue, loaded_xc_search, loaded_yc_search, loaded_zc_search, loaded_creation_time, loaded_particle_mass, loaded_current_time = numpy.loadtxt("{0}/particle_{1}_location.txt".format(filesys_location, withParticleIDValue), unpack=True)
				for value in range(len(loaded_file_number)):
					#pulled_file_number = loaded_file_number[value]
					#pulled_file_number = str(pulled_file_number)[:-2]
					pulled_file_number = file_number(str(loaded_file_number[value])[:-2])
					if file_number == pulled_file_number:
						print file_number, pulled_file_number
						print loaded_ParticleIDValue[0], withParticleIDValue
						xc_search = loaded_xc_search[value]
						yc_search = loaded_yc_search[value]
						zc_search = loaded_zc_search[value]
						#this_creation_time = loaded_creation_time[value]
						print 'xc = ', xc_search, 'yc = ', yc_search,'zc = ', zc_search
						break
					else:
						continue
				print 'Have found the search location, now searching for max density.'

			elif First_non_particle and not withRestart:
				# We haven't restarted the run and are now at the first file that doesn't have particles.
				xc_search = xc
				yc_search = yc
				zc_search = zc
				First_non_particle = False
			else:
				# i.e. we're continuing to search backwards already
				xc_search = xc_new
				yc_search = yc_new
				zc_search = zc_new
			Search_radius_Sphere = 0.5

			file = prefix+"{0:05d}/info_{0:05d}.txt".format(i)
			print file
			sinkfile = prefix+"{0:05d}/sink_{0:05d}.info".format(i)
			print sinkfile
			pf = yt.load(file)

			#pf = load("{0}{1:04d}".format( plt_prefix, i))
			#sp = pf.h.sphere((xc_search, yc_search, zc_search), Search_radius_Sphere/pf['pc'])
			#dens = sp["Density"]
			#max= dd.quantities["MaxLocation"]("Density")
			#max = sp.quantities["MaxLocation"]("Density")
			#maxDens = max[0]/3e-22
			#maxLoc = numpy.array(max[2:5])/3e18
			# set the new density peak location
			#xc_new = max[2]
			#yc_new = max[3]
			#zc_new = max[4]
			current_time = pf.current_time
			particleMass = 0.0
			this_particle_id = int(withParticleIDValue)
			dd = pf.all_data()
			position_search = YTArray( [xc_search, yc_search, zc_search], "cm")
			sp = pf.sphere(position_search, (Search_radius_Sphere, 'pc'))
			max_location = sp.quantities.max_location("Density")
			# this convert from YTArray to just almost normal numbers
			xc_new = max_location[2].in_cgs()#.to_ndarray()
			yc_new = max_location[3].in_cgs()#.to_ndarray()
			zc_new = max_location[4].in_cgs()#.to_ndarray()

			vxc = 0
			vyc = 0
			vzc = 0
			print 'These are the center of the search coordinates: ', xc_search, yc_search, zc_search
			print 'These are the Max Density coordinates: ',xc_new, yc_new, zc_new
			print withParticleIDValue, this_creation_time
			print 'Passing to the analysis script'
			pass_to_analysis_script(pf, xc_new, yc_new, zc_new, this_particle_id, this_creation_time, current_time, radiusMin, radiusSphere, particleMass, file_number, sinkfile, file)

#############################################################
#########_______MAIN_TRACING_BACKWARDS_IN_TIME_END______#####
#############################################################
#
###########################################################
########_______MAIN_TRACING_FORWARDS_IN_TIME_______########
###########################################################
def Main_tracing_forwards():
	for i in range(args.start,args.end,args.step) :	
		# First, check that the two files exist:
		#FLASH
		if withFLASH4:
			if not (flash_file_existance_chk):
				continue
			print 'On File: {0}{1:04d}'.format(plt_prefix, i)
		else:
			print 'On File: output_{0:05d}'.format(i)
		# Create the 5 digit file number so that we can reference it for plotting
		file_number = standardize_file_number(i)
		#print four_digit_file_number
		file = prefix+"{0:05d}/info_{0:05d}.txt".format(i)
		print file
		#file_exist = glob.glob(file)
		sinkfile = prefix+"{0:05d}/sink_{0:05d}.info".format(i)
		print sinkfile
		pf = yt.load(file)
		#particle_ID, partMass, xstar, ystar, zstar = Obtain_particles(sinkfile)
		particle_ID, partMass, xstar, ystar, zstar, vxstar, vystar, vzstar = Obtain_particles(sinkfile)
		# Need creation time, and current time from RAMSES
		print particle_ID
		print particle_ID.size

		if particle_ID.size == 0:
			print 'Particles have not been created yet'
			#print 'Setting the particle list to include the specified particle id.'
			#particle_ID = [withParticleIDValue]
			print 'Exiting.'
			sys.exit()
		elif particle_ID.size == 1:
			if particle_ID == withParticleIDValue or (withAllParticles):
				xc = xstar	 
				yc = ystar	 
				zc = zstar
				print xc, yc, zc, partMass
				xc = pf.quan(xc, "cm")
				yc = pf.quan(yc, "cm")
				zc = pf.quan(zc, "cm")
				print xc, yc, zc, partMass
				particleMass = partMass
				this_particle_id = int(particle_ID)
				print this_particle_id
						#this_creation_time = creation_time[j]
				this_creation_time = 0.0
				particle_age = 0.0
				current_time = 0.0
				if (withShell) or (withSmallSphere) or (withShellSphere):
					print 'going to python script'
					pass_to_analysis_script(pf, xc, yc, zc, this_particle_id, this_creation_time, current_time, radiusMin, radiusSphere, particleMass, file_number, sinkfile, file)
					#fileout="/home/m/murray/dwmurray/scratch/test_ramses_mhd2/{1}_{2:04d}_{3}_{4}.out".format(filesys_location, out_prefix, i, compare_file, this_particle_id)
					#copyfile(sinkfile, '/home/m/murray/dwmurray/scratch/ramses_jeans_{0}/sink_{1:05d}.info'.format(filesys_location, i))
					#copyfile(file, '/home/m/murray/dwmurray/scratch/ramses_jeans_{0}/info_{1:05d}.txt'.format(filesys_location, i))
					#print fileout
					#getRadialProfile_py(pf, xc, yc, zc, this_particle_id, this_creation_time, current_time, fileout, radiusMin, radiusSphere, particleMass)

		elif particle_ID.size >= 1:
			half_point = xstar.size/2
			creation_time = 0.0
			current_time = 0.0
			
			if (withFirstHalf):
				print 'will only do the first half of the list of particles.'
				for j in range(xstar.size/2 + 1):
					part_lst_num = j
					print 'On particle:', part_lst_num + 1, 'of:', xstar.size, ' in File: {0}{1:04d}'.format(plt_prefix, i)
					print particle_ID[part_lst_num], withParticleIDValue
					if particle_ID[part_lst_num] == withParticleIDValue or (withAllParticles):
						xc, yc, zc, vxc, vyc, vzc, this_particle_id, this_creation_time, particleMass = RAMSES_obtain_individual_part_attributes(part_lst_num, file_number, pf, particle_ID, xstar, ystar, zstar, vxstar, vystar, vzstar, partMass, creation_time, current_time)
					else:
						print 'This particle ID does not match, ' + \
						    'may not be in this octant.'
						print 'Check the particle_location.txt file for ' + \
						    'a list of particles that are in this quadrant.'
						continue
					print 'Passing to the analysis script'
					pass_to_analysis_script(pf, xc, yc, zc, this_particle_id, this_creation_time, current_time, radiusMin, radiusSphere, particleMass, file_number, sinkfile, file)

			elif (withSecondHalf):
				print 'will only do the Second half of the list of particles.'
				for j in range(xstar.size/2 + 1):
					part_lst_num = j+half_point
					print 'On particle:', part_lst_num + 1, 'of:', xstar.size, ' in File: {0}{1:04d}'.format(plt_prefix, i)
					if part_lst_num + 1 > xstar.size:
						# This would query past the end length of xstar.size.
						break
					print particle_ID[part_lst_num], withParticleIDValue
					if particle_ID[part_lst_num] == withParticleIDValue or (withAllParticles):
						xc, yc, zc, vxc, vyc, vzc, this_particle_id, this_creation_time, particleMass = RAMSES_obtain_individual_part_attributes(part_lst_num, file_number, pf, particle_ID, xstar, ystar, zstar, vxstar, vystar, vzstar, partMass, creation_time, current_time)
#						xc, yc, zc, vxc, vyc, vzc, this_particle_id, this_creation_time, particleMass = obtain_individual_part_attributes(part_lst_num, file_number, pf, particle_ID, xp, yp, zp, vxp, vyp, vzp, partMass, creation_time, current_time)
					else:
						print 'This particle ID does not match, ' + \
						    'may not be in this octant.'
						print 'Check the particle_location.txt file for ' + \
						    'a list of particles that are in this quadrant.'
						continue
					print 'Passing to the analysis script'
					pass_to_analysis_script(pf, xc, yc, zc, this_particle_id, this_creation_time, current_time, radiusMin, radiusSphere, particleMass, file_number, sinkfile, file)
			elif (withParallel) : 
				print 'run in parallel.'
				parallel_size = int(math.ceil(1.0 * xstar.size/size)) 
				for j in range(rank, parallel_size, size):
					part_lst_num = j
					print 'On particle:', part_lst_num + 1, 'of:', xstar.size, ' in File: {0}{1:04d}'.format(plt_prefix, i)
					if part_lst_num  >= xstar.size:
						# This would query past the end length of xstar.size.
						break
					print particle_ID[part_lst_num], withParticleIDValue
					if particle_ID[part_lst_num] == withParticleIDValue or (withAllParticles):
						xc, yc, zc, vxc, vyc, vzc, this_particle_id, this_creation_time, particleMass = RAMSES_obtain_individual_part_attributes(part_lst_num, file_number, pf, particle_ID, xstar, ystar, zstar, vxstar, vystar, vzstar, partMass, creation_time, current_time)
#						xc, yc, zc, vxc, vyc, vzc, this_particle_id, this_creation_time, particleMass = obtain_individual_part_attributes(part_lst_num, file_number, pf, particle_ID, xp, yp, zp, vxp, vyp, vzp, partMass, creation_time, current_time)
					else:
						print 'This particle ID does not match, ' + \
						    'may not be in this octant.'
						print 'Check the particle_location.txt file for ' + \
						    'a list of particles that are in this quadrant.'
						continue
					print 'Passing to the analysis script'
					pass_to_analysis_script(pf, xc, yc, zc, this_particle_id, this_creation_time, current_time, radiusMin, radiusSphere, particleMass, file_number, sinkfile, file)
			else:
				print 'will do the whole list of particles.'
				for j in range(xstar.size):
					part_lst_num = j
					print 'On particle:', part_lst_num + 1, 'of:', xstar.size, ' in File: {0}{1:04d}'.format(plt_prefix, i)
					print particle_ID[part_lst_num], withParticleIDValue
					if particle_ID[part_lst_num] == withParticleIDValue or (withAllParticles):
						xc, yc, zc, vxc, vyc, vzc, this_particle_id, this_creation_time, particleMass = RAMSES_obtain_individual_part_attributes(part_lst_num, file_number, pf, particle_ID, xstar, ystar, zstar, vxstar, vystar, vzstar, partMass, creation_time, current_time)
						#xc, yc, zc, vxc, vyc, vzc, this_particle_id, this_creation_time, particleMass = obtain_individual_part_attributes(part_lst_num, file_number, pf, particle_ID, xp, yp, zp, vxp, vyp, vzp, partMass, creation_time, current_time)
					else:
						print 'This particle ID does not match, ' + \
						    'may not be in this octant.'
						print 'Check the particle_location.txt file for ' + \
						    'a list of particles that are in this quadrant.'
						continue
					print 'Passing to the analysis script'
					pass_to_analysis_script(pf, xc, yc, zc, this_particle_id, this_creation_time, current_time, radiusMin, radiusSphere, particleMass, file_number, sinkfile, file)
				
#
#
#			for j in range(xstar.size) :
#				print particle_ID[j], withParticleIDValue
#				if particle_ID[j] == withParticleIDValue or (withAllParticles):
#					xc = xstar[j]	 
#					yc = ystar[j]	 
#					zc = zstar[j]
#					print xc, yc, zc, partMass[j]
#					xc = pf.quan(xc, "cm")
#					yc = pf.quan(yc, "cm")
#					zc = pf.quan(zc, "cm")
#				#vxc = vxstar[j]
#				#vyc = vystar[j]
#				#vzc = vzstar[j]
#					print xc, yc, zc, partMass[j]
#					particleMass = partMass[j]
#					this_particle_id = particle_ID[j]
#					print this_particle_id
#				#this_creation_time = creation_time[j]
#					this_creation_time = 0.0
#					particle_age = 0.0
#					current_time = 0.0
#				#particle_age = current_time - this_creation_time
#				#particle_age = particle_age / (numpy.pi*1e7)
#					print 'On particle:', j + 1, 'of:', xstar.size, 'in File: {0}{1:05d}'.format(prefix, i)
#				else:
#					print 'This particle ID does not match, ' + \
#					    'may not be in this octant.'
#					print 'Check the particle_location.txt file for ' + \
#					    'a list of particles that are in this quadrant.'
#					continue
#
#
#		with open("particle_location.txt", "a") as particle_location:
#			particle_location.write("Particle ID's in this file: " + str(particle_ID))
#			particle_location.write('\n')
#			particle_location.close()
#

###########################################################
########____END_MAIN_TRACING_FORWARDS_IN_TIME_______#######
###########################################################


def Main_tracing_forwards_no_particles():
	for i in range(args.start,args.end,args.step) :	
		# First, check that the two files exist:
		if not (flash_file_existance_chk):
			continue
		print 'On File: {0}{1:04d}'.format(plt_prefix, i)
		# Create the 4 digit file number so that we can reference it for plotting
		file_number = standardize_file_number(i)
		#print four_digit_file_number
		# Use this for before star particle formation.
		# Make sure getRadialProfile_yt has bulk set to sphere and not particle Velocity
		pf = yt.load("{0}{1:04d}".format( plt_prefix, i))
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
		current_time = pf.current_time
		#this_particle_id = int(withParticleIDValue)
		this_particle_id = int(42)
		this_creation_time = 0.0
		particle_age = 0.0
		this_particle_Mass = 0.0
		print 'Passing to the analysis script'
		pass_to_analysis_script(pf, xc, yc, zc, this_particle_id, this_creation_time, current_time, radiusMin, radiusSphere, particleMass, file_number, sinkfile, file)

###########################################################
###########################################################

#pf = load("BB_hdf5_plt_cnt_0096")

# compressive case

import argparse
parser = argparse.ArgumentParser(description = "start number to end number, stride length, reduction method")

parser.add_argument('start', metavar='N1', type=int, help ='Start value for hdf5 files')
parser.add_argument('end', metavar='N2', type=int, help='End Value for hdf5 files, note runs until End-1')
parser.add_argument('step', metavar='N3', type=int, help='Stride length')
parser.add_argument('ParticleID', metavar='N4', type=int, nargs='?', default=42, help='Particle ID you want to reduce.')
parser.add_argument('filesyslocation_to_output2', metavar='N5', type=str, nargs='?', default='filsys', help='file_sys location')
parser.add_argument('--smallsphere', action='store_true')
parser.add_argument('--bigsphere', action='store_true')
parser.add_argument('--particle', action='store_true')
parser.add_argument('--noparticle', action='store_true')
parser.add_argument('--shell', action='store_true')
parser.add_argument('--shellsphere', action='store_true')
parser.add_argument('--chk', action='store_true')
parser.add_argument('--allparticles', action='store_true')
parser.add_argument('--backwards', action='store_true')
parser.add_argument('--restart', action='store_true')
parser.add_argument('--first', action='store_true')
parser.add_argument('--second', action='store_true')
parser.add_argument('--parallel', action='store_true')
parser.add_argument('--FLASH4', action='store_true')
#parser.add_argument('--partrestart', action='store_true')
args = parser.parse_args()

# Some of These settings eventually go to yt calls
withSmallSphere = args.smallsphere
withBigSphere = args.bigsphere
withParticle = args.particle
withNoParticle = args.noparticle
withCheckpoint = args.chk
withParticleIDValue = args.ParticleID
withAllParticles = args.allparticles
withBackTracking = args.backwards
withRestart = args.restart
withFilesyslocation = args.filesyslocation_to_output2
withFirstHalf = args.first
withSecondHalf = args.second
withParallel = args.parallel
withFLASH4 = args.FLASH4
# These calls eventually go to the python written outputs
withShell = args.shell
withShellSphere = args.shellsphere

Sphere_Bulk = False
Particle_Bulk = False
NoParticle_Bulk = False

Shell_Bulk = False
ShellSphere_Bulk = False

if (withParallel):
	from mpi4py import MPI

	comm = MPI.COMM_WORLD
	size = comm.Get_size()
	rank = comm.Get_rank()

#filesys_location = str(withFilesyslocation)
#filesys_location = '/home/m/murray/dwmurray/scratch/test-ramses/{0}/python_output'.format(withFilesyslocation)
filesys_location = '/Users/dwmurray/Work/Ramses_hydro/{0}'.format(withFilesyslocation)
if (withSmallSphere):
	compare_file = 'smallsphere'
	Sphere_Bulk = True
	radiusMin=1e-3
	radiusSphere=0.15

if (withBigSphere):
	compare_file = 'bigsphere'
	radiusMin=1e-3
	radiusSphere=3.0
	Sphere_Bulk = True
	Bulk_sphere_radius = 3.0

if (withParticle):
	compare_file = 'part'
	radiusMin=1e-3
	radiusSphere=3.0
	Particle_Bulk = True

if (withShell):
	compare_file = 'shell'
	radiusMin=1e-3
	radiusSphere=3.0
	Shell_Bulk = True


if (withShellSphere):
	compare_file = 'shellsphere'
	radiusMin=1e-3
	radiusSphere=3.0
	ShellSphere_Bulk = True

# These are universal for FLASH
if (withCheckpoint):
	plt_prefix = "BB_hdf5_chk_"
	part_prefix = "BB_hdf5_part_"
else:
	plt_prefix = "BB_hdf5_plt_cnt_"
	part_prefix = "BB_hdf5_part_"

# This list is for any I matricies that are singular,
# or have issues with numpy.linalg.inv
# If they work with the penrose psuedo invert, I'd like to know
Penrose_matrix_particles = []
prefix = "output_"
out_prefix = "rad_profile"
quad = os.getcwd()[-5:]

# This is what we store in the txt files
# (zip_store_file_number, zip_store_particle_ID, zip_store_xc, zip_store_yc, zip_store_zc, zip_store_creation_time, zip_store_particle_mass, zip_store_current_time)
#with open("particle_location.txt", "a") as particle_location:
#	particle_location.write('Frame #' + " ")
#	particle_location.write('Particle ID' + " ")
#	particle_location.write('X value' + " ")
#	particle_location.write('Y value' + " ")
#	particle_location.write('Z value' + " ")
#	particle_location.write('Creation Time (yr)' + " ")
#	particle_location.write('Particle Mass (g)' + " ")
#	particle_location.write('Current Time (yr)' + " ")
#	particle_location.write('\n')
#	particle_location.close()
#
# The file looks like this:
#'{out_prefix}{framestep}_{compare_file}_{particle_number}.out'                                                                                               
# i.e. rad_profile_0218_part_000.out   
# compare file options area: bigsphere, smallsphere, part, nopart, shell, shellsphere
if (withBackTracking):
	trackbackwards()
elif (withNoParticle):
	Main_tracing_forwards_no_particles()
else:
	Main_tracing_forwards()

print 'These particles needed a pseudo-inversion, not sure if you should trust them'
print Penrose_matrix_particles
