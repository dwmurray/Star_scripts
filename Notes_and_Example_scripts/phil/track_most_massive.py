import matplotlib
matplotlib.use("Agg")

from yt.mods import * 
import numpy 
import argparse

bins = 100
meanDens = 3e-22
i = 0
Msun = 1.99e33
G = 6.67e-8
parsec = 3.09e18

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

		
def getRadialProfile(pf, xc, yc, zc, fileout="rad_profile.out", radiusSphere=4.0, particleMass = 0.0) : 
	sp = pf.h.sphere([xc,yc,zc], radiusSphere/pf['pc'])

	sp2 = pf.h.sphere([xc,yc,zc], 0.5/pf['pc'])
	com = sp2.quantities["CenterOfMass"]()
	xcom = com[0]
	ycom = com[1]
	zcom = com[2]


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
	numpy.savetxt(fileout, zip(rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, particleMass*numpy.ones(mTbin.size)), fmt="%15.9E")

parser = argparse.ArgumentParser(description = "start")
parser.add_argument('Start', metavar='N1', type=int)
parser.add_argument('End', metavar='N2', type=int)
parser.add_argument('Last', metavar='N3', type=int)
args = parser.parse_args()
startIndex = args.Start
endIndex = args.End
lastIndex = args.Last

# first get the masses of the most massive particles
prefix = "BB_hdf5_plt_cnt_"
partTrack = 5
pf = load(prefix + "{0:04d}".format(lastIndex))

data = pf.h.all_data()
partInd  = data["particle_index"]
partMass = data["ParticleMassMsun"]
partX = data["particle_position_x"]
partY = data["particle_position_y"]
partZ = data["particle_position_z"]

sortInd = numpy.argsort(partMass)[::-1] # order by most massive first

partMass = partMass[sortInd][0:partTrack]
partInd = partInd[sortInd][0:partTrack]
partX = partX[sortInd][0:partTrack]
partY = partY[sortInd][0:partTrack]
partZ = partZ[sortInd][0:partTrack]

numpy.savetxt( "most_massive_stars.data", zip( partInd, partMass, partX, partY, partZ))
MassHistory = numpy.zeros( [endIndex+1 - startIndex,partTrack])
TimeHistory = numpy.zeros( endIndex+1-startIndex)
for fileIndex in range(startIndex,endIndex+1) : 
	pf = load("{0}{1:04d}".format(prefix,fileIndex))
	data = pf.h.all_data()
	pindex = data["particle_index"]
	pmass = data["ParticleMassMsun"]
	px = data["particle_position_x"]
	py = data["particle_position_y"]
	pz = data["particle_position_z"]

	TimeHistory[ fileIndex -startIndex ] = fileIndex
	for j in range(partTrack) :  # 10 most massive particles
		for i in range(pindex.size) : 
			if ( partInd[j] == pindex[i]) :
				MassHistory[fileIndex-startIndex,j] =  pmass[i]
				xc = px[i]
				yc = py[i]
				zc = pz[i]
				getRadialProfile(pf, xc, yc, zc, fileout="radial_profile_{0:02d}_{1:03d}.data".format(j, fileIndex),radiusSphere=3.0, particleMass = pmass[i])
				
				

numpy.savetxt( "most_massive.data", zip( TimeHistory, MassHistory[:,0], MassHistory[:,1], MassHistory[:,2], MassHistory[:,3], MassHistory[:,4]))
#, MassHistory[:,5], MassHistory[:,6], MassHistory[:,7],  MassHistory[:,8], MassHistory[:,9]))

	#TrackAnalysis(pf, partInd, fileoutPrefix = "Track_{0:04d}".format(fileIndex))  
 
