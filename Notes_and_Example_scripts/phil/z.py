import matplotlib
matplotlib.use("Agg")

from yt.mods import * 
import numpy 
import argparse
import shutil
import math 

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

		getRadialProfile( pf, xp, yp, zp, fileout="{0}_{1:03d}.data".format(fileoutPrefix,i), particleMass=mps[i])

		
def getRadialProfile(pf, xc, yc, zc, fileout="rad_profile.out", radiusSphere=4.0, particleMass = 0.0) : 
	sp = pf.h.sphere([xc,yc,zc], radiusSphere/pf['pc'])

	sp2 = pf.h.sphere([xc,yc,zc], 0.01/pf['pc'])
	bulk = sp.quantities["BulkVelocity"]()
	sp.set_field_parameter("bulk_velocity", bulk)
	rp = BinnedProfile1D( sp, 50, "Radiuspc", 1e-3, radiusSphere, log_space=True)
	rp.add_fields("Density", weight="CellVolume")
	maxDens = rp["Density"][0]
	rp.add_fields("RadialVelocity")  
	rp.add_fields("VelocityMagnitude")
	rp.add_fields("SpecificAngularMomentumX")
	rp.add_fields("SpecificAngularMomentumY")
	rp.add_fields("SpecificAngularMomentumZ")
	rp.add_fields("CellMassMsun", accumulation=True,weight=None)
	mT = rp["CellMassMsun"][1:]+particleMass
        rho = rp["Density"][1:]
	radius = rp["Radiuspc"][1:]
	vk = numpy.sqrt(G*mT*Msun/(parsec*radius))
	vk = vk[1:]
	vr = rp["RadialVelocity"][1:]
	lx = rp["SpecificAngularMomentumX"][1:]
	ly = rp["SpecificAngularMomentumY"][1:]
	lz = rp["SpecificAngularMomentumZ"][1:]
	vmag = rp["VelocityMagnitude"][1:]
        vrms = numpy.sqrt(vmag**2. - vr**2.)
	l = numpy.sqrt( lx*lx + ly*ly + lz*lz)
	numpy.savetxt(fileout, zip(radius, vr, vrms, vk, vmag, mT, rho, l), fmt="%15.9E")


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
parser.add_argument('Last', metavar='N3', type=int)
args = parser.parse_args()
startIndex = args.Start
endIndex = args.End
lastIndex = args.Last

# first get the masses of the most massive particles
prefix = "BB_hdf5_plt_cnt_"
partTrack =100 
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

dr0 = 0.1*parsec
for i in range(partX.size) :
	position = [partX[i], partY[i], partZ[i]]
	dr = numpy.sqrt((partX - position[0])**2 + (partY - position[1])**2 + (partZ - position[2])**2)
	mass = partMass[ dr < dr0].sum()
	print partMass[i], mass, numpy.array(position)/3.08e18
	#makeDiskFrame(prefix, lastIndex, position[0], position[1], position[2], radiusSphere = 3.0,
	#	fileout = "Disk_{0:04d}_{1:04d}.png".format(lastIndex, i))  
	#makeDiskFrame(prefix, lastIndex, position[0], position[1], position[2], radiusSphere = 0.04,
	#	fileout = "movie_disk_{0:04d}_{1:04d}.png".format(lastIndex, i))  
	AccretionProfile(pf, position[0], position[1], position[2], fileout="Accretion_{0:03d}_0.5pc.npz".format(i), radiusSphere=0.5)
	AccretionProfile(pf, position[0], position[1], position[2], fileout="Accretion_{0:03d}_0.02pc.npz".format(i), radiusSphere=0.02, deltaR=0.002)

for fileIndex in range(startIndex,endIndex+1) : 
	pf = load("{0}{1:04d}".format(prefix,fileIndex))
	TrackAnalysis(pf, partInd, fileoutPrefix = "Track_{0:04d}".format(fileIndex))  
 
