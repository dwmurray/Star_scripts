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

def TrackAnalysis(pf, partInd, fileoutPrefix = "track") : 
	data = pf.h.all_data()

	partIndices = data["particle_index"]
	xps = data["particle_position_x"]
	yps = data["particle_position_y"]
	zps = data["particle_position_z"]
	mps = data["ParticleMassMsun"]

	for i in range(partInd.size) :
		pIndex = partInd[i]

		if not pIndex in partIndices :
			break		
		
		index = numpy.where(partIndices == partInd[i])[0][0]

		xp = xps[index]
		yp = yps[index]
		zp = zps[index]
		mp = mps[index]

		getRadialProfile( pf, xp, yp, zp, fileout="{0}_{1:03d}.data".format(fileoutPrefix,i))

def getVelocityPowerSpec( pf, xc, yc, zc, cubeSize = 1.0,sigma=1.0,fileout="power_spec.pdf") :
	global window
	global windowInitialized
        level = 5
        left_edge = numpy.array([xc - 0.5*cubeSize/pf['pc'], yc - 0.5*cubeSize/pf['pc'], zc - 0.5*cubeSize/pf['pc']])
        right_edge = numpy.array([xc + 0.5*cubeSize/pf['pc'], yc + 0.5*cubeSize/pf['pc'], zc + 0.5*cubeSize/pf['pc']])

        dims = (16*2**level * (right_edge-left_edge)/(pf.domain_right_edge-pf.domain_left_edge)).astype(int)
        cube = pf.h.covering_grid(level=level,left_edge=left_edge,right_edge=right_edge,dims=dims)
        #cube = pf.h.periodic_region_strict(left_edge=left_edge,right_edge=right_edge)
        startime = time.time()
        print 'start vel ' + str(time.time()-startime)
        vel = cube["VelocityMagnitude"]

        size = dims[0]
        step = 2.0/size
        arr = numpy.arange(-1.0+step*0.5, 1.0, step)
        arr2= arr*arr
        sigma2inv = cubeSize*cubeSize/(sigma*sigma)
	if( windowInitialized == False) : 
	        window = numpy.zeros(dims)
		windowInitialized = True
        	for i in range(dims[0]) :
                	for j in range(dims[1]) :
                        	for k in range(dims[2]) :
                                	window[i,j,k] = math.exp(-(arr2[i]+arr2[j]+arr2[k])*sigma2inv)

	vel = vel*window
        print 'start filter ' + str(time.time()-startime)
        ffdat = ff.FourierFilter(vel)

        print 'end filter ' + str(time.time()-startime)
        bins = ffdat.bins[:ffdat.nx]

        print 'start power ' + str(time.time()-startime)
        power = na.abs(na.fft.fftn(vel)/vel.size)**2
        print 'end power' + str(time.time()-startime)
        spec = na.zeros(ffdat.nx)
        for i in range(ffdat.nx):
                for j in range(ffdat.nx) :
                        for k in range(ffdat.nx) :
                                bin = int(math.floor(math.sqrt(i*i + j*j + k*k)))
                                if(bin < ffdat.nx) :
                                        spec[bin] = spec[bin] + power[i,j,k]

        norm = (spec*bins*bins).sum()
        print 'end spec ' + str(time.time()-startime)
        numpy.savetxt(fileout, zip(bins, spec))


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
	density = sp["Density"]

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
	nbin = numpy.zeros(bins)
	vrbin = numpy.zeros(bins)
	vxbin = numpy.zeros(bins)
	vybin = numpy.zeros(bins)
	vzbin = numpy.zeros(bins)

	# get center of mass for each shell
#	for i in range(r.size) : 
#		index = int(r[i]*bins/(radiusSphere*parsec))
#		if(index < bins) :
#			mbin[index] = mbin[index] + cellMass[i]
#			nbin[index] = nbin[index] + 1
#			vxbin[index] = vxbin[index] + vx[i]*cellMass
#			vybin[index] = vybin[index] + vy[i]*cellMass
#			vzbin[index] = vzbin[index] + vz[i]*cellMass
#	vxbin = vxbin/mbin
#	vybin = vybin/mbin
#	vzbin = vzbin/mbin

	mbin = numpy.zeros(bins)
	nbin = numpy.zeros(bins)
	# get the radial velocity
	for i in range(r.size) : 
		index = int(r[i]*bins/(radiusSphere*parsec))
		if(index < bins) :
			mbin[index] = mbin[index] + cellMass[i]
			nbin[index] = nbin[index] + 1
			vr = ((vx[i]-vxbin[index])*x[i] + (vy[i]-vybin[index])*y[i] + (vz[i]-vzbin[index])*z[i])/r[i]*cellMass[i]
			vrbin[index] = vrbin[index] + vr
		
	vrbin = vrbin/mbin
	#vrbin = vrbin/nbin

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
			densbin[index] = densbin[index] + density[i]

	vrmsbin = numpy.sqrt(vrmsbin/mbin)
	vrmsnbin = numpy.sqrt(vrmsnbin/nbin)
	vmagbin = numpy.sqrt(vmagbin/mbin)
	vmagnbin = numpy.sqrt(vmagnbin/nbin)
	densbin = densbin/nbin

	# get the Kepler velocity
	mTbin = mbin
	mTbin[0] = mTbin[0] + particleMass
	for i in range(1,mTbin.size) : 
		mTbin[i] = mTbin[i] + mTbin[i-1]
		
	rbin = radiusSphere*numpy.arange(bins)/bins
	vKbin = numpy.sqrt( G*mTbin*Msun/((rbin + 0.01)*parsec))
	numpy.savetxt(fileout, zip(rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin,densbin), fmt="%15.9E")

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

				index = int(bins*dist/(radiusSphere*parsec))

				if(index < bins) :
					vrmsbin[index] = vrmsbin[index] + vrms
					nbin[index] = nbin[index] + 1
		

	vrmsbin = numpy.sqrt(vrmsbin/nbin)
	rbin = radiusSphere*numpy.arange(bins)/bins

	numpy.savetxt(fileout, zip(rbin, vrmsbin), fmt="%12.7e")


def mslice( pf, fileout = 'junk.png',field="Density", center="max", width=4) :
	newImage = Image.new("RGB", (1070,2000))
	print fileout 
	p = SlicePlot(pf, "x", field, center=center,width=(width,"pc") )
	p.set_zlim("Density", 1e-23,1e-19)
	p.annotate_velocity(factor=16)
	
	pid = os.getpid()
        p.save("junk{0:06d}".format(pid))

	p = SlicePlot(pf, "y", field, center=center, width=(4,"pc") )
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

def AccretionProfile(pf, xc, yc, zc, fileout="rad_profile.out", radiusSphere=4.0, deltaR = 0.05, particleMass = 0.0) : 
	sp = pf.h.sphere([xc,yc,zc], (radiusSphere+deltaR)/pf['pc'])
	numPhi = 30
	numTheta = 15

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

	rho = sp["Density"]

	r = numpy.sqrt(x*x + y*y + z*z)

	vx = sp["x-velocity"]
	vy = sp["y-velocity"]
	vz = sp["z-velocity"]

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

	pi = math.pi

	rhodotSph = numpy.zeros([numTheta, numPhi])
	rhoSph = numpy.zeros([numTheta, numPhi])
	iSph = numpy.zeros([numTheta, numPhi])
	aSph = numpy.zeros([numTheta, numPhi])

	dTheta = pi/numTheta*0.5
	dPhi = 2*pi/numTheta*0.5
	aSph[0,:] = (1.-math.cos(dTheta))*dPhi
	for i in range(1,numTheta) : 
		theta = pi * i/numTheta + dTheta
		theta1 = theta - dTheta
		theta2 = theta + dTheta
		aSph[i,:] = -(math.cos(theta2) - math.cos(theta1))*dPhi

	# get the radial velocity
	for i in range(r.size) : 
		if(abs(r[i] - radiusSphere*parsec)/parsec < deltaR) : 
			theta = math.acos( z[i]/r[i])
			rxy = math.sqrt(x[i]**2 + y[i]**2)
			phi = math.atan2(y[i]/rxy, x[i]/rxy)
			if(rxy < 0.1) : 
				phi = 0.
			rhodot = (vx[i]*x[i] + vy[i]*y[i] + vz[i]*z[i])*rho[i]/r[i]
			indexPhi = int((phi+pi)*numPhi/(2.*pi))
			indexTheta = int((theta)*numTheta/(pi))
			rhodotSph[indexTheta,indexPhi] = rhodotSph[indexTheta,indexPhi] + rhodot
			rhoSph[indexTheta,indexPhi] = rhoSph[indexTheta,indexPhi] + rho[i]
			iSph[indexTheta,indexPhi] = iSph[indexTheta,indexPhi] + 1
	
	rhodotSph = rhodotSph/numpy.maximum(iSph,numpy.ones([numTheta,numPhi]))
	rhoSph = rhoSph/numpy.maximum(iSph,numpy.ones([numTheta,numPhi]))
	mdotSph = rhodotSph * aSph
	numpy.savez( fileout, rhodotSph = rhodotSph, rhoSph=rhoSph, aSph=aSph)

parser = argparse.ArgumentParser(description = "start")
parser.add_argument('Start', metavar='N1', type=int)
parser.add_argument('End', metavar='N2', type=int)
args = parser.parse_args()
startIndex = args.Start
endIndex = args.End

# first get the masses of the most massive particles
prefix = "BB_hdf5_plt_cnt_"

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

for fileIndex in range(endIndex,startIndex,-1) : 
	pf = load("{0}{1:04d}".format(prefix,fileIndex))

	partIndices = None
	xps = None
	yps = None
	zps = None
	mps = None
	
	# get all the star particles first 
	if("particle_index" in pf.h.derived_field_list) : 
		all_data = pf.h.all_data()
	
		partIndices = all_data["particle_index"]
		xps = all_data["particle_position_x"]
		yps = all_data["particle_position_y"]
		zps = all_data["particle_position_z"]
		mps = all_data["ParticleMassMsun"]

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
					if (r < 0.5*3e18) : 
						found = True
				if( not found) : 
					xfound.append(xp)
					yfound.append(yp)
					zfound.append(zp)
					indfound.append(index)
					print "Found density at ({0:5.2f}, {1:5.2f}, {2:5.2f})".format(xp/3e18, yp/3e18, zp/3e18)
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
		maxLoc = numpy.array(max[2:5])/3e18
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
		#pc.plots[-1].data.write_out("track_density_{0:03d}_{1:02d}.data".format(i,-int(ncenters[i])), format='%15.9E', bin_style='left')
		#pc.save()

#		mslice( pf, fileout = 'frame{0:02d}_{1:02d}.png'.format(i,int(ncenters[i])),field="Density", center=[xc,yc,zc], width=4) 
		#AccretionProfile(pf, xc, yc, zc, fileout="Accretion_{0:03d}_t={1:02d}.npz".format(i,int(ncenters[i])), radiusSphere=0.5)

	#	getVelocityPowerSpec(pf, xc, yc, zc, cubeSize=8.0, sigma=2.0, fileout="track_power_spec{0:03d}_t={1:02d}.data".format(i,int(ncenters[i])))
		#getRadialProfile(pf, xc, yc, zc, fileout="TrackDensity_{0:03d}_t={1:02d}.data".format(i,-int(ncenters[i])), radiusSphere=3.0)
	#	getLarsonsLaw(pf, xc, yc, zc, fileout="track_larson_{0:03d}_t={1:02d}.data".format(i,int(ncenters[i])), radiusSphere=4.0, randPoints=False)
		#shutil.move("{0}{1:04d}_Profile1D_0_Radiuspc_Density.png".format(prefix, fileIndex), "track_hi_density_{0:04d}.png".format(fileIndex))
	

for i in range(indcenters.size) : 
	print "{0} {1}".format(i,int(indcenters[i]))
 
pf = load("{0}{1:04d}".format(prefix,84))
for i in range(0,10): 
	xc = random.random()*4.8e19
	yc = random.random()*4.8e19
	zc = random.random()*4.8e19
	#getLarsonsLaw(pf, xc, yc, zc, fileout="track_larson_random_{0:03d}.data".format(i), radiusSphere=4.0, randPoints=False)
#	getVelocityPowerSpec(pf, xc, yc, zc, cubeSize=16.0, fileout="track_power_spec_random_point{0:04d}.data".format(i))
	#getVelocityPowerSpec(pf, xc, yc, zc, cubeSize=8.0, sigma=2.0, fileout="track_power_spec_random_point{0:04d}.data".format(i))
