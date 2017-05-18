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

bins = 50
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

				index = int(bins*dist/(radiusSphere*parsec))

				if(index < bins) :
					vrmsbin[index] = vrmsbin[index] + vrms
					nbin[index] = nbin[index] + 1
		

	vrmsbin = numpy.sqrt(vrmsbin/nbin)
	rbin = radiusSphere*numpy.arange(bins)/bins

	numpy.savetxt(fileout, zip(rbin, vrmsbin), fmt="%12.7e")

parser = argparse.ArgumentParser(description = "start")
parser.add_argument('Start', metavar='N1', type=int)
parser.add_argument('End', metavar='N2', type=int)
args = parser.parse_args()
startIndex = args.Start
endIndex = args.End

# first get the masses of the most massive particles
prefix = "BB_hdf5_plt_cnt_"

initialized = False

for fileIndex in range(endIndex,startIndex,-1) : 
	pf = load("{0}{1:04d}".format(prefix,fileIndex))
	if( initialized ) : 
		data = pf.h.sphere([xc,yc,zc], 0.5/pf['pc'])
	else : 
		data = pf.h.all_data()
		initialized = True 
	max= data.quantities["MaxLocation"]("Density")
	maxDens = max[0]/3e-22
	maxLoc = numpy.array(max[2:5])/3e18
	print fileIndex, maxDens, maxLoc
	xc = max[2]
	yc = max[3]
	zc = max[4]
	#pc = PlotCollection(pf, center = [xc,yc,zc])
	#pc.add_profile_sphere( 5.0, "pc", ["Radiuspc", "Density"], weight="CellMassMsun",x_bins=100,x_log=False)
	#pc.plots[-1].data.write_out("track_density_{0:04d}.data".format(fileIndex), format='%15.9E', bin_style='left')
	#pc.save()
	getVelocityPowerSpec(pf, xc, yc, zc, cubeSize=4.0, sigma=0.5, fileout="track_power_spec{0:04d}_sigma0.5.data".format(fileIndex))
	#getRadialProfile(pf, xc, yc, zc, fileout="TrackDensity_{0:03d}.data".format(fileIndex), radiusSphere=3.0)
	#getLarsonsLaw(pf, xc, yc, zc, fileout="track_larson_{0:03d}.data".format(fileIndex), radiusSphere=3.0, randPoints=False)
	#shutil.move("{0}{1:04d}_Profile1D_0_Radiuspc_Density.png".format(prefix, fileIndex), "track_hi_density_{0:04d}.png".format(fileIndex))
	


pf = load("{0}{1:04d}".format(prefix,120))
for i in range(0,5): 
	xc = random.random()*4.8e19
	yc = random.random()*4.8e19
	zc = random.random()*4.8e19
	#getLarsonsLaw(pf, xc, yc, zc, fileout="track_larson_random_{0:03d}.data".format(i), radiusSphere=3.0, randPoints=False)
	#getVelocityPowerSpec(pf, xc, yc, zc, cubeSize=16.0, fileout="track_power_spec_random_point{0:04d}.data".format(i))
	getVelocityPowerSpec(pf, xc, yc, zc, cubeSize=4.0, sigma=0.5, fileout="track_power_spec_random_point{0:04d}_sigma0.5.data".format(i))
