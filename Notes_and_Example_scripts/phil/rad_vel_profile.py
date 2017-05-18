import matplotlib
matplotlib.use("Agg")

from yt.mods import *
import matplotlib.pyplot as p
import numpy

pf = load("BB_hdf5_plt_cnt_0134")
dd = pf.h.all_data()
dd["particle_posx"]
xparr = dd["particle_posx"]
yparr = dd["particle_posy"]
zparr = dd["particle_posz"]
partMassArr = dd["ParticleMassMsun"]

radiusSphere = 4.0
meanDens = 3e-22
i = 0
Msun = 1.99e33
G = 6.67e-8
parsec = 3.09e18
for xp,yp,zp,partMass in zip(xparr,yparr,zparr,partMassArr) :
	sp = pf.h.sphere([xp,yp,zp],radiusSphere/pf['pc'])
	bulk = sp.quantities["BulkVelocity"]()
	sp.set_field_parameter("bulk_velocity", bulk)
	rp = BinnedProfile1D( sp, 50, "Radiuspc", 0., radiusSphere, log_space=False)
	rp.add_fields("Density")
	maxDens = rp["Density"][0]
	if( maxDens > 1000.*meanDens) :
		rp.add_fields("RadialVelocity")  
		rp.add_fields("VelocityMagnitude")
		rp.add_fields("CellMassMsun", accumulation=True,weight=None)
		vk = numpy.sqrt(G*(rp["CellMassMsun"]+partMass)*Msun/(parsec*rp["Radiuspc"]))
		i = i + 1
		radius = rp["Radiuspc"][1:]
		vk = vk[1:]
		p.plot(radius, -rp["RadialVelocity"][1:]/vk, label="$-v_r/v_K$ $M_*={0:4.2f}$".format(partMass))
	        vrms = numpy.sqrt(rp["VelocityMagnitude"]**2. - rp["RadialVelocity"]**2.)
		p.plot(radius, vrms[1:]/vk, linestyle="--", label="$v_{\\rm rms}/v_K$")
		#p.plot(rp["Radiuspc"], rp["CellMassMsun"])

		p.xlabel( "$r$ [pc]")
		p.ylabel( "$v/v_K$")
		p.legend(loc='best')
		p.savefig("radial_velocity_{0:02d}.pdf".format(i))
		p.clf()
