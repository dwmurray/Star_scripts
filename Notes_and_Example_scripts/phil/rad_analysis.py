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
logRhoMin = -4
logRhoMax = 4.

def getRadialProfile(pf, xc, yc, zc, fileout="rad_profile.out", radiusSphere=3.0, particleMass = 0.0) : 
	sp = pf.h.sphere([xc,yc,zc], radiusSphere/pf['pc'])

	x = sp["x"] - xc
	y = sp["y"] - yc
	z = sp["z"] - zc
	r = numpy.sqrt(x*x + y*y + z*z)
	cellMass = sp["CellMassMsun"]
	dens = sp["Density"]
	cellVolume = sp["CellVolume"]

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
	rhobin = numpy.zeros(bins)
	volbin = numpy.zeros(bins)

	# get the radial velocity
	for i in range(r.size) : 
		index = int(r[i]*bins/(radiusSphere*parsec))
		if(index < bins) :
			mbin[index] = mbin[index] + cellMass[i]
			volbin[index] = volbin[index] + cellVolume[i]
			vr = (vx[i]*x[i] + vy[i]*y[i] + vz[i]*z[i])*cellMass[i]/r[i]
			vrbin[index] = vrbin[index] + vr
			rhobin[index] = rhobin[index] + dens[i]
		
	vrbin = vrbin/mbin
	rhobin = rhobin/volbin

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
	numpy.savetxt(fileout, zip(rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, rhobin), fmt="%15.9E")

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



#pf = load("BB_hdf5_plt_cnt_0096")

# compressive case
pf = load("BB_hdf5_plt_cnt_0077")
dd = pf.h.all_data()
dd.quantities[["particle_posx"]
xparr = dd["particle_posx"]
yparr = dd["particle_posy"]
zparr = dd["particle_posz"]
partMassArr = dd["ParticleMassMsun"]

radiusSphere = 3.0
xc = 8.*parsec
yc = 8.*parsec
zc = 8.*parsec

# for full case
#getDensityExclude( pf, xc,yc,zc, xparr, yparr, zparr, fileoutrho_PDF="rho_PDFExclude_3pc.out", radiusExclude=3.0, exclude=True) 
#getDensityFunctions(pf, 0., 0., 0., fileoutrho_r="rho_r0.out", fileoutrho_PDF="rho_PDF_Full.out".format(i), 
#	radiusSphere=5.0, allData = True) 

# for compressive case
getDensityFunctions(pf, 0., 0., 0., fileoutrho_r="rho_r0.out", fileoutrho_PDF="rho_PDF_Compressive.out".format(i), 
	radiusSphere=5.0, allData = True) 
