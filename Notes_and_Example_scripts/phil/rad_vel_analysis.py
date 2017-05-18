import matplotlib
matplotlib.use("Agg")

from yt.mods import *
import matplotlib.pyplot as p

import math
import numpy
import random

bins = 50
meanDens = 3e-22
i = 0
Msun = 1.99e33
G = 6.67e-8
parsec = 3.09e18

def getRadialProfile(pf, xc, yc, zc, fileout="rad_profile.out", radiusSphere=3.0, particleMass = 0.0) : 
	sp = pf.h.sphere([xc,yc,zc], radiusSphere/pf['pc'])

	x = sp["x"] - xc
	y = sp["y"] - yc
	z = sp["z"] - zc
	r = numpy.sqrt(x*x + y*y + z*z)
	cellMass = sp["CellMassMsun"]

	vx = sp["x-velocity"]
	vy = sp["y-velocity"]
	vz = sp["z-velocity"]

	# First figure out the bulk velocity
	Mtot = cellMass.sum()
	vbx = (vx*cellMass).sum()/Mtot
	vby = (vy*cellMass).sum()/Mtot
	vbz = (vz*cellMass).sum()/Mtot

	# correct for bulk velocity
	vx = vx - vbx
	vy = vy - vby
	vz = vz - vbz

	# now get the radial velocity

	# first get the central velocity
	sp2 = pf.h.sphere([xc,yc,zc], 0.05/pf['pc'])
	vcx = sp2["x-velocity"]
	vcy = sp2["y-velocity"]
	vcz = sp2["z-velocity"]
	cMass = sp2["CellMassMsun"]
	cMassTot = cMass.sum()

	vcx = (cMass*vcx).sum()/cMassTot - vbx
	vcy = (cMass*vcy).sum()/cMassTot - vby
	vcz = (cMass*vcz).sum()/cMassTot - vbz
	
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

	# get the rms velocity
	for i in range(r.size) : 
		index = int(r[i]*bins/(radiusSphere*parsec))
		if(index < bins) :
			mbin[index] = mbin[index] + cellMass[i]
			nbin[index] = nbin[index] + 1
			rad = r[i]
			if(rad < 0.01) :
				 rad = 0.01
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



pf = load("BB_hdf5_plt_cnt_0144")
dd = pf.h.all_data()
dd["particle_posx"]
xparr = dd["particle_posx"]
yparr = dd["particle_posy"]
zparr = dd["particle_posz"]
partMassArr = dd["ParticleMassMsun"]

radiusSphere = 3.0

for i in range(partMassArr.size) : 
	if (partMassArr[i] > 30.) : 
		xc = xparr[i]
		yc = yparr[i]
		zc = zparr[i]

		getRadialProfile(pf, xc, yc, zc, fileout="rad_profile_{0:02d}.out".format(i), 
			radiusSphere=radiusSphere,particleMass=partMassArr[i])
		getLarsonsLaw(pf, xc, yc, zc, fileout="larson_{0:02d}.out".format(i), 
			radiusSphere=3.0,trials=3,randPoints = True) 

