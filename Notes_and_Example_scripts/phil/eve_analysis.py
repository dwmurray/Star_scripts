from yt.mods import *

#for i in range(len(index)):
#for i in np.arange(1, len(index)-1):
#i = 0
#if i != -2:

pf = load("BB_hdf5_plt_cnt_0088")
data = pf.h.all_data()

partIndices = data["particle_index"]
xps = data["particle_position_x"]
yps = data["particle_position_y"]
zps = data["particle_position_z"]
mps = data["ParticleMassMsun"]

sp = pf.h.sphere(p_pos[i,:], 5./pf["pc"])
p = pf.h.find_particles_by_type(4)

if type(p) != types.NoneType:
    particles = True
else:
    particles = False

#com = sp.quantities["CenterOfMass"](use_particles=particles)
if particles:
    max = sp.quantities["MaxLocation"]("Matter_Density")
else:
    max = sp.quantities["MaxLocation"]("Density")
com = max[2:5] #maximum density cell
rad = np.linspace(0., 5., 51)

#calculate central velocity
if particles:
    masslabel = "TotalMassMsun"
else:
    masslabel = "CellMassMsun"
spc = pf.h.sphere(com, 0.05/pf["pc"])
vcx = spc.quantities["WeightedAverageQuantity"]("x-velocity",masslabel)
vcy = spc.quantities["WeightedAverageQuantity"]("y-velocity",masslabel)
vcz = spc.quantities["WeightedAverageQuantity"]("z-velocity",masslabel)
spc = 0. #save memory

##initialize velocity arrays    
    #Infall Velocity
vr = np.zeros(len(rad)-1)
vr_m = np.zeros(len(rad)-1)

    #Vmag
vmag = np.zeros(len(rad)-1) #vol-weighted
vmag_m = np.zeros(len(rad)-1) #mass-weighted

#free-fall
vff = np.zeros(len(rad)-1)

#load sphere
sp = pf.h.sphere(com, 5./pf["pc"])

#radius
sp_r = sp["Radiuspc"]
x = sp["x"]-com[0]
y = sp["y"]-com[1]
z = sp["z"]-com[2]

#bulk velocity
vbx, vby, vbz = sp.get_field_parameter("bulk_velocity")

#velocity
vx = sp["x-velocity"]-vcx
vy = sp["y-velocity"]-vcy
vz = sp["z-velocity"]-vcz

#mass
mass = sp[masslabel]

for j in range(1, len(rad)):
    index = np.where((sp_r >= rad[j-1])*(sp_r < rad[j]) > 0)[0]

    vr[j-1] = np.mean(sp["RadialVelocity"][index])

    if vr[j-1] > 0: vr[j-1] = 0
    vrx = vr[j-1]*(sp["x"][index]-com[0])/(rad[j]/pf["pc"])
    vry = vr[j-1]*(sp["y"][index]-com[1])/(rad[j]/pf["pc"])
    vrz = vr[j-1]*(sp["z"][index]-com[2])/(rad[j]/pf["pc"])

    vmag[j-1] = np.mean(np.sqrt((vx[index]-vrx)**2.+(vy[index]-vry)**2.+(vz[index]-vrz)**2.))
    vmag_m[j-1] = np.sum(mass[index]*np.sqrt((vx[index]-vrx)**2.+(vy[index]-vry)**2.+(vz[index]-vrz)**2.))/np.sum(mass[index])

    vff[j-1] = np.sqrt(consts.G*mass[np.where(sp_r < rad[j])[0]].sum()*consts.Msun/rad[j]/consts.pc2cm)

#take only infall velocity
vr = np.where(vr < 0, np.abs(vr), 0.)
np.savetxt(outdir+"starN%05i_DD%04i_vprof.txt" %(stari, picki), zip(vr, vmag, vmag_m, vff, rad))

