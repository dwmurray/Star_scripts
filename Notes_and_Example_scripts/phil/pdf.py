import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as pl
from yt.mods import *
import numpy 

mean_dens = 3e-22
pf = load("BB_hdf5_plt_cnt_0122")
region = pf.h.all_data()

pc = PlotCollection(pf)
pc.add_profile_object( region, ["Density", "CellMassMsun"], weight=None)

dens_mass = pc.plots[-1].data["Density"]/mean_dens
mass =  pc.plots[-1].data["CellMassMsun"]

pc.add_profile_object( region, ["Density", "CellVolume"], weight=None)
dens_volume = pc.plots[-1].data["Density"]/mean_dens
volume =  pc.plots[-1].data["CellVolume"]

mass = mass/mass.sum()
volume = volume/volume.sum()

pl.plot( dens_mass, mass, 'b-', label="mass")
pl.plot( dens_volume, volume, 'b--', label="volume")
pl.xscale('log')
pl.yscale('log')
pl.ylim(1e-8, 1e-1)
pl.legend(loc='best')
pl.xlabel("$\\rho/\\rho_0$")
pl.ylabel("PDF")
pl.savefig("prob.pdf")
