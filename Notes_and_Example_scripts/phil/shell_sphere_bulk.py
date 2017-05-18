import matplotlib
matplotlib.use("Agg")

from yt.mods import *
import matplotlib.pyplot as p

import math
import numpy
import random

bins = 35
meanDens = 3e-22
i = 0
Msun = 1.99e33
G = 6.67e-8
parsec = 3.09e18
logRhoMin = -4
logRhoMax = 4.

def getRadialProfile(pf, xc, yc, zc, fileout="rad_profile.out", radiusMin=1e-3, radiusSphere=3.0, particleMass = 0.0) : 
	sp = pf.h.sphere((xc, yc, zc), radiusSphere/pf['pc'])

	x = sp["x"] - xc
	y = sp["y"] - yc
	z = sp["z"] - zc
	r = numpy.sqrt(x*x + y*y + z*z)
	# this is in Solar Masses
	cellMass = sp["CellMassMsun"]
	# This is in cgs
	dens = sp["Density"]
	# This is in cm**3
	cellVolume = sp["CellVolume"]

	# These are in cm/s
	vx = sp["x-velocity"]
	vy = sp["y-velocity"]
	vz = sp["z-velocity"]

	# This is the total Mass
	Mtot = cellMass.sum()

	lx = y*vz - vy*z
	ly = z*vx - x*vz
	lz = x*vy - y*vx

	mbin = numpy.zeros(bins)
	vxbin = numpy.zeros(bins)
	vybin = numpy.zeros(bins)
	vzbin = numpy.zeros(bins)

	vx_bulk_bin = numpy.zeros(bins)
	vy_bulk_bin = numpy.zeros(bins)
	vz_bulk_bin = numpy.zeros(bins)

	vx_sphere_bin = numpy.zeros(bins)
	vy_sphere_bin = numpy.zeros(bins)
	vz_sphere_bin = numpy.zeros(bins)

	rhobin = numpy.zeros(bins)
	volbin = numpy.zeros(bins)
	vrbin = numpy.zeros(bins)
	mdotbin = numpy.zeros(bins)
	menc = numpy.zeros(bins)
#####################################

	# Find the Bulk Velocity of each shell, prior to removing from each cell in the shell
	print "Calculating the Bulk Velocity"
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
			vx_bulk_bin[index] = vx_bulk_bin[index] + (vx[i]*cellMass[i])
			vy_bulk_bin[index] = vy_bulk_bin[index] + (vy[i]*cellMass[i])
			vz_bulk_bin[index] = vz_bulk_bin[index] + (vz[i]*cellMass[i])

        for shell_test in range(1,bins):
		#menc[shell_test] = menc[shell_test] + menc[shell_test-1]
		menc[shell_test] = mbin[shell_test] + menc[shell_test-1]
                vx_bulk_bin[shell_test] = vx_bulk_bin[shell_test] + vx_bulk_bin[shell_test-1] 
                vy_bulk_bin[shell_test] = vy_bulk_bin[shell_test] + vy_bulk_bin[shell_test-1] 
                vz_bulk_bin[shell_test] = vz_bulk_bin[shell_test] + vz_bulk_bin[shell_test-1] 
		#vz_bulk_bin[shell_test] = vz_bulk_bin[shell_test]/mbin[shell_test]
                #for shell_sphere in range(shell_test):
                #        vx_sphere_bin[shell_test] = vx_bulk_bin[shell_sphere]
                #        vy_sphere_bin[shell_test] = vy_sphere_bin[shell_test] + vy_bulk_bin[shell_sphere]
                #        vz_sphere_bin[shell_test] = vz_sphere_bin[shell_test] + vz_bulk_bin[shell_sphere]
		#print vx_sphere_bin[shell_test]
		
	vx_sphere_bin = vx_bulk_bin/menc
	vy_sphere_bin = vy_bulk_bin/menc
	vz_sphere_bin = vz_bulk_bin/menc

	
#####################################
	#mbin = numpy.zeros(bins)
	print "Obtaining vr"
	# get the radial velocity
	lgradiusMin = math.log10( radiusMin)
	lgradiusSph = math.log10( radiusSphere)

	for i in range(r.size) : 
		if( r[i]/parsec < radiusMin) :
			continue
		index = int((math.log10(r[i]/parsec)-lgradiusMin)*bins/(lgradiusSph - lgradiusMin))
		if(index >= 0 and index < bins) :
			# and the Volume of each shell
			volbin[index] = volbin[index] + cellVolume[i]

			# This vr is actually Mdot * distance
			vr_x_mass = ((vx[i] - vx_sphere_bin[index])*x[i] + (vy[i] - vy_sphere_bin[index])*y[i] + (vz[i] - vz_sphere_bin[index])*z[i])*cellMass[i]/r[i]

			rhobin[index] = rhobin[index] + dens[i]*cellVolume[i]
			#mdotbin[index] = mdotbin[index] + vr_x_mass/r[i] # vr is mdot*r right now

			vrbin[index] = vrbin[index] + vr_x_mass


	vrbin = vrbin/mbin
	print vrbin
	#print mbin
	#rhobin = mbin/volbin
	rhobin = rhobin/volbin

	mbin = numpy.zeros(bins)
	vrmsbin = numpy.zeros(bins)
	vrmsnbin = numpy.zeros(bins)
	vmagbin = numpy.zeros(bins)
	vmagnbin = numpy.zeros(bins)
	angXbin = numpy.zeros(bins)
	angYbin = numpy.zeros(bins)
	angZbin = numpy.zeros(bins)
	nbin = numpy.zeros(bins)

	print "Obtaining vrms"
	# get the rms velocity
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
			vrmsn = (vx[i]-vrbin[index]*x[i]/rad)**2. + (vy[i]-vrbin[index]*y[i]/rad)**2. + (vz[i]-vrbin[index]*z[i]/rad)**2.
			vrms = vrmsn*cellMass[i]
			vmagn = (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i])
			vmag = vmagn*cellMass[i]
			vrmsbin[index] = vrmsbin[index] + vrms
			vrmsnbin[index] = vrmsnbin[index] + vrmsn
			vmagbin[index] = vmagbin[index] + vmag
			vmagnbin[index] = vmagnbin[index] + vmagn
			angXbin[index] = lx[i] * cellMass[i] + angXbin[index]
			angYbin[index] = ly[i] * cellMass[i] + angYbin[index]
			angZbin[index] = lz[i] * cellMass[i] + angZbin[index]

	vrmsbin = numpy.sqrt(vrmsbin/mbin)
	vrmsnbin = numpy.sqrt(vrmsnbin/nbin)
	vmagbin = numpy.sqrt(vmagbin/mbin)
	vmagnbin = numpy.sqrt(vmagnbin/nbin)

	# get the Kepler velocity
	print "Calculating the Kepler Velocity"
	mTbin = mbin
	mTbin[0] = mTbin[0] + particleMass
	for i in range(1,mTbin.size) : 
		mTbin[i] = mTbin[i] + mTbin[i-1]
		angXbin[i] = angXbin[i] + angXbin[i-1]
		angYbin[i] = angYbin[i] + angYbin[i-1]
		angZbin[i] = angZbin[i] + angZbin[i-1]
		
	angXbin = angXbin/mTbin
	angYbin = angYbin/mTbin
	angZbin = angZbin/mTbin

	norm = numpy.sqrt(angXbin**2 + angYbin**2 + angZbin**2)
	angXbin = angXbin/norm
	angYbin = angYbin/norm
	angZbin = angZbin/norm
	lgrbin = lgradiusMin + (lgradiusSph-lgradiusMin)*numpy.arange(bins)/bins
	rbin = 1e1**lgrbin
	vKbin = numpy.sqrt( G*mTbin*Msun/((rbin)*parsec))
	numpy.savetxt(fileout, zip(rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, rhobin, mdotbin, norm, angXbin, angYbin, angZbin), fmt="%15.9E")

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

import argparse

#plt_prefix = "BB_hdf5_chk_"
plt_prefix = "BB_hdf5_plt_cnt_"
out_prefix = "sphere_bulk_test"
parser = argparse.ArgumentParser(description = "start number to end number")

parser.add_argument('start', metavar='N1', type=int)
parser.add_argument('end', metavar='N2', type=int)
parser.add_argument('step', metavar='N2', type=int)

args = parser.parse_args()

for i in range(args.start,args.end,args.step) :	

	pf = load("{0}{1:04d}".format( plt_prefix, i))

	dd = pf.h.all_data()
	xp = dd["particle_posx"]
	yp = dd["particle_posy"]
	zp = dd["particle_posz"]
	partMass = dd["ParticleMassMsun"]

	#max= data.quantities["MaxLocation"]("Density")
	#maxDens = max[0]/3e-22
	#maxLoc = numpy.array(max[2:5])/3e18
	#xc = max[2]
	#yc = max[3]
	#zc = max[4]
	
	for j in range(xp.size) :
		xc = xp[j]	 
		yc = yp[j]	 
		zc = zp[j]	 
		getRadialProfile(pf,xc,yc,zc, fileout="{0}{1:04d}_part_{2:03d}.out".format( out_prefix, i, j), radiusSphere=3.0, particleMass = partMass[j])
