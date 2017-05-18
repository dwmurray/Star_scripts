import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as p

from yt.mods import *
import numpy

# first define a moment of inertia derived quantity

def _MomentOfInertia(data, xc=0., yc=0., zc=0.) : 
	x = data["x"] - xc
	y = data["y"] - yc
	z = data["z"] - zc
	m = data["CellMassMsun"]
	Ixx = (m*y*y + m*z*z).sum()
	Iyy = (m*x*x + m*z*z).sum()
	Izz = (m*x*x + m*y*y).sum()
	Ixy = -(m*x*y).sum()
	Ixz = -(m*x*z).sum()
	Iyz = -(m*y*z).sum()

	return Ixx, Iyy, Izz, Ixy, Ixz, Iyz

def _combMomentOfInertia(data, Ixx, Iyy, Izz, Ixy, Ixz, Iyz) : 
	IxxTot = Ixx.sum()
	IyyTot = Iyy.sum()
	IzzTot = Izz.sum()
	IxyTot = Ixy.sum()
	IxzTot = Ixz.sum()
	IyzTot = Iyz.sum()

	return numpy.array([[IxxTot, IxyTot, IxzTot],
		           [IxyTot,IyyTot, IyzTot],
                           [IxzTot, IyzTot, IzzTot]])

add_quantity("MomentOfInertia", function=_MomentOfInertia,
             combine_function=_combMomentOfInertia, n_ret = 6)

fileIndex = 93
# first get the masses of the most massive particles
prefix = "BB_hdf5_plt_cnt_"
partStart = 0
partTrack = 1
pf = load(prefix + "{0:04d}".format(fileIndex))


data = pf.h.all_data()
partInd  = data["particle_index"]
partMass = data["ParticleMassMsun"]
partX = data["particle_position_x"]
partY = data["particle_position_y"]
partZ = data["particle_position_z"]

print partMass
sortInd = numpy.argsort(partMass)[::-1] # order by most massive first

partMass = partMass[sortInd][partStart:partTrack]
partInd = partInd[sortInd][partStart:partTrack]
partX = partX[sortInd][partStart:partTrack]
partY = partY[sortInd][partStart:partTrack]
partZ = partZ[sortInd][partStart:partTrack]


for x,y,z,m in zip(partX,partY,partZ,partMass) :
	rbin = []
	I1 = []
	I2 = []
	Mt = []
	for r in numpy.arange(0.25, 5., 0.25) : 
		sp = pf.h.sphere([x,y,z], radius=r/pf['pc'])
		com = sp.quantities["CenterOfMass"]()
		I = sp.quantities["MomentOfInertia"](com[0], com[1], com[2])
		#I = sp.quantities["MomentOfInertia"](x, y, z)
		eigValues, eigVectors = numpy.linalg.eig(I)
		rbin.append(r)
		sortEig = numpy.sort(eigValues)/eigValues.max()
		I1.append(sortEig[0])
		I2.append(sortEig[1])
		Mt.append(sp.quantities["TotalMass"]())
	
		print r, eigValues/eigValues.max(), Mt[-1]
	I1 = numpy.array(I1)
	I2 = numpy.array(I2)
	Mt = numpy.array(Mt)
	rbin = numpy.array(rbin)

	p.loglog(rbin, I1, label="$I_1$")
	p.loglog(rbin, I2, label="$I_2$")
#	p.loglog(rbin, Mt/Mt.max(), label="$M(r)$")
#	p.loglog(numpy.arange(0.3,1.0,0.1),0.1*numpy.arange(0.3,1.0,0.1))
#	p.loglog(numpy.arange(2.0,5.0,0.1),0.1*numpy.arange(2.0,5.0,0.1)**2.0/4.)
	p.legend(loc=4)
	p.xlabel("r")
	#p.ylabel("$I_1/I_3,\ I_2/I_3,\ M(r)/M(r=5\,{\\rm pc})$")
	p.ylabel("$I_1/I_3,\ I_2/I_3$")
	p.xlim(0.25,5.0)
	p.ylim(0.1, 1.0)
	p.savefig("geometry.pdf")
