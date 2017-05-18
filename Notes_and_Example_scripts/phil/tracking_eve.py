import matplotlib
matplotlib.use("Agg")
import sys
import math

sys.path.append('./scripts/joishi/fourier_tools/fourier_tools/')

import fourier_filter as ff

from yt.mods import * 
import numpy 
import argparse
import shutil
import matplotlib.pyplot as pl
import random
import Image
import os

bins = 80
meanDens = 3e-22
i = 0
Msun = 1.99e33
G = 6.67e-8
parsec = 3.09e18
window = None
windowInitialized = False


def getRadialProfile(pf, xc, yc, zc, fileout="rad_profile.out", radiusSphere=4.0, particleMass = 0.0) : 
	sp = pf.h.sphere([xc,yc,zc], radiusSphere/pf['pc'])

	sp2 = pf.h.sphere([xc,yc,zc], 0.05/pf['pc'])
	com = sp2.quantities["CenterOfMass"]()
	xcom = com[0]
	ycom = com[1]
	zcom = com[2]

#	xc = xcom
#	yc = ycom
#	zc = zcom

	x = sp["x"] - xcom
	y = sp["y"] - ycom
	z = sp["z"] - zcom

	r = numpy.sqrt(x*x + y*y + z*z)

	vx = sp["x-velocity"]
	vy = sp["y-velocity"]
	vz = sp["z-velocity"]

	cellMass = sp["CellMassMsun"]

	# correct for bulk velocity
	vx = vx
	vy = vy
	vz = vz

	# now get the radial velocity

	# first get the central velocity
	sp2 = pf.h.sphere([xc,yc,zc], 0.05/pf['pc'])
	vcx = sp2["x-velocity"]
	vcy = sp2["y-velocity"]
	vcz = sp2["z-velocity"]
	cMass = sp2["CellMassMsun"]
	cMassTot = cMass.sum()

	vcx = (cMass*vcx).sum()/cMassTot
	vcy = (cMass*vcy).sum()/cMassTot
	vcz = (cMass*vcz).sum()/cMassTot
	
	# get velocities relative to center
	vx = vx - vcx
	vy = vy - vcy
	vz = vz - vcz

	mbin = numpy.zeros(bins)
	vrbin = numpy.zeros(bins)

	# get the radial velocity
	for i in range(r.size) : 
		index = int(r[i]*bins/(radiusSphere*parsec))
		if(index < bins) :
			mbin[index] = mbin[index] + cellMass[i]
			vr = (vx[i]*x[i] + vy[i]*y[i] + vz[i]*z[i])*cellMass[i]/r[i]
			vrbin[index] = vrbin[index] + vr
		
	vrbin = vrbin/mbin

	mbin = numpy.zeros(bins)
	vrmsbin = numpy.zeros(bins)
	vrmsnbin = numpy.zeros(bins)
	vmagbin = numpy.zeros(bins)
	vmagnbin = numpy.zeros(bins)
	nbin = numpy.zeros(bins)
	densbin = numpy.zeros(bins)

	# get the rms velocity
	for i in range(r.size) : 
		index = int(r[i]*bins/(radiusSphere*parsec))
		if(index < bins) :
			mbin[index] = mbin[index] + cellMass[i]
			nbin[index] = nbin[index] + 1
			rad = r[i]
			if(rad < 0.01/16.0) : #change for ENZO
				 rad = 0.01/16.0 
			vrmsn = (vx[i]-vrbin[index]*x[i]/rad)**2. + (vy[i]-vrbin[index]*y[i]/rad)**2. + (vz[i]-vrbin[index]*z[i]/rad)**2.
			vrms = vrmsn*cellMass[i]
			vmagn = (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i])
			vmag = vmagn*cellMass[i]
			vrmsbin[index] = vrmsbin[index] + vrms
			vrmsnbin[index] = vrmsnbin[index] + vrmsn
			vmagbin[index] = vmagbin[index] + vmag
			vmagnbin[index] = vmagnbin[index] + vmagn

	vrmsbin = numpy.sqrt(vrmsbin/mbin)
	vrmsnbin = numpy.sqrt(vrmsnbin/nbin)
	vmagbin = numpy.sqrt(vmagbin/mbin)
	vmagnbin = numpy.sqrt(vmagnbin/nbin)

	# get the Kepler velocity
	mTbin = mbin
	mTbin[0] = mTbin[0] + particleMass
	for i in range(1,mTbin.size) : 
		mTbin[i] = mTbin[i] + mTbin[i-1]
		
	rbin = radiusSphere*numpy.arange(bins)/bins
	vKbin = numpy.sqrt( G*mTbin*Msun/((rbin + 0.01)*parsec))
	numpy.savetxt(fileout, zip(rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin), fmt="%15.9E")


parser = argparse.ArgumentParser(description = "start")
parser.add_argument('Start', metavar='N1', type=int)
parser.add_argument('End', metavar='N2', type=int)
args = parser.parse_args()
startIndex = args.Start
endIndex = args.End

filename = "DD{0:04d}/data{0:04d}"

initialized = False

# list of centers 
xcenters = numpy.array([])
ycenters = numpy.array([])
zcenters = numpy.array([])
ncenters = numpy.array([])
indcenters = numpy.array([])

starIndices = None
xstars = None
ystars = None
zstars = None
rfound = 0.5*3e18 # for FLASH
rfound = 0.5/16.0 # for enzo (code units)

for fileIndex in range(endIndex,startIndex,-1) : 
	pf = load(filename.format(fileIndex))

	partIndices = None
	xps = None
	yps = None
	zps = None
	mps = None
	
	# get all the star particles first 
	if("particle_index" in pf.h.derived_field_list) : 
		all_data = pf.h.all_data()
		
		validParticles = all_data["particle_type"]==4	
		partIndices = all_data["particle_index"][validParticles]
		xps = all_data["particle_position_x"][validParticles]
		yps = all_data["particle_position_y"][validParticles]
		zps = all_data["particle_position_z"][validParticles]
		mps = all_data["ParticleMassMsun"][validParticles]

	if( initialized and partIndices != None) :
		# find missing star particles
		xfound = []
		yfound = []
		zfound = []
		indfound = []
	
		print "looking for star particles"
		for i in range(starIndices.size) :
			index = starIndices[i] 
			if not index in partIndices : 	
				# check if something is close enough
				xp = xstars[i]
				yp = ystars[i]
				zp = zstars[i]

				found = False
				for partIndex,x,y,z in zip(partIndices,xps,yps,zps) : 
					r = math.sqrt( (x-xp)*(x-xp) + (y-yp)*(y-yp) + (z-zp)*(z-zp))
					if (r < rfound) : 
						found = True
				if( not found) : 
					xfound.append(xp)
					yfound.append(yp)
					zfound.append(zp)
					indfound.append(index)
					print "Found density at ({0:5.2f}, {1:5.2f}, {2:5.2f})".format(xp, yp, zp)
		if( xfound) : # found something now append to list
			xcenters = numpy.append(xcenters,xfound)
			ycenters = numpy.append(ycenters,yfound)
			zcenters = numpy.append(zcenters,zfound)
			ncenters = numpy.append(ncenters,numpy.zeros(len(xfound)))
			indcenters = numpy.append(indcenters, indfound)
		# redefine star indices
		starIndices = partIndices
		xstars = xps
		ystars = yps
		zstars = zps

	elif( partIndices != None) : 
		starIndices = partIndices
		xstars = xps
		ystars = yps
		zstars = zps
		xcenters = numpy.array([])
		ycenters = numpy.array([])
		zcenters = numpy.array([])
		ncenters = numpy.array([])
		initialized = True 

	print ncenters
	print xcenters 
	for i in range(ncenters.size) :  
		xc = xcenters[i] 
		yc = ycenters[i] 
		zc = zcenters[i] 
		data = pf.h.sphere([xc,yc,zc], 0.5/pf['pc'])
		max= data.quantities["MaxLocation"]("Density")
		maxDens = max[0]/3e-22
		maxLoc = numpy.array(max[2:5])
		print fileIndex, maxDens, maxLoc
		xc = max[2]
		yc = max[3]
		zc = max[4]
		
		# redefine the new xcenters
		xcenters[i] = xc
		ycenters[i] = yc
		zcenters[i] = zc
		ncenters[i] = ncenters[i] + 1

		print "Working on {0} at position ({1:5.2f}, {2:5.2f}, {3:5.2f})".format(i,xc/3e18, yc/3e18, zc/3e18)
	
		#pc = PlotCollection(pf, center = [xc,yc,zc])
		#pc.add_profile_sphere( 5.0, "pc", ["Radiuspc", "Density"], weight="CellVolume",x_bins=150,x_log=False)
		#pc.plots[-1].data.write_out("track_density_{0:03d}_{1:02d}.data".format(i,int(ncenters[i])), format='%15.9E', bin_style='left')
		#pc.save()

		#mslice( pf, fileout = 'frame{0:02d}_{1:02d}.png'.format(i,int(ncenters[i])),field="Density", center=[xc,yc,zc], width=4) 
		#AccretionProfile(pf, xc, yc, zc, fileout="Accretion_{0:03d}_t={1:02d}.npz".format(i,int(ncenters[i])), radiusSphere=0.5)
		#getVelocityPowerSpec(pf, xc, yc, zc, cubeSize=8.0, sigma=2.0, fileout="track_power_spec{0:04d}.data".format(fileIndex))
		getRadialProfile(pf, xc, yc, zc, fileout="EveLeeDensity_{0:03d}_t={1:02d}.data".format(i,int(ncenters[i])), radiusSphere=3.0)
	#	getLarsonsLaw(pf, xc, yc, zc, fileout="track_larson_{0:03d}_t={1:02d}.data".format(i,int(ncenters[i])), radiusSphere=4.0, randPoints=False)
		#shutil.move("{0}{1:04d}_Profile1D_0_Radiuspc_Density.png".format(prefix, fileIndex), "track_hi_density_{0:04d}.png".format(fileIndex))
	

for i in range(indcenters.size) : 
	print "{0} {1}".format(i,int(indcenters[i]))
 
#pf = load("{0}{1:04d}".format(prefix,90))
for i in range(0,10): 
	xc = random.random()*4.8e19
	yc = random.random()*4.8e19
	zc = random.random()*4.8e19
	#getLarsonsLaw(pf, xc, yc, zc, fileout="track_larson_random_{0:03d}.data".format(i), radiusSphere=4.0, randPoints=False)
	#getVelocityPowerSpec(pf, xc, yc, zc, cubeSize=16.0, fileout="track_power_spec_random_point{0:04d}.data".format(i))
	#getVelocityPowerSpec(pf, xc, yc, zc, cubeSize=8.0, sigma=2.0, fileout="track_power_spec_random_point{0:04d}.data".format(i))
