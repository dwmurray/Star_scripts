import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as pl
from yt.mods import *
import numpy 

def plotPDF( pf, posx, posy, posz, rmax=5.0, label="5 pc", color='blue') : 
	
	mean_dens = 3e-22
	region = pf.h.sphere([posx,posy,posz], rmax/pf['pc'])

	pc = PlotCollection(pf)
	pc.add_profile_object( region, ["Density", "CellMassMsun"], weight=None)

	dens_mass = pc.plots[-1].data["Density"]/mean_dens
	mass =  pc.plots[-1].data["CellMassMsun"]

	p = pc.add_profile_object( region, ["Density", "CellVolume"], weight=None)

	dens_volume = pc.plots[-1].data["Density"]/mean_dens
	volume =  pc.plots[-1].data["CellVolume"]

	mass = mass/mass.sum()
	volume = volume/volume.sum()

	pl.plot( dens_mass, mass, '-', color=color, label="{0} mass".format(label))
	pl.plot( dens_volume, volume, '--', color=color, label="{0} volume".format(label))

	lgDens = numpy.log10(dens_volume)
	lgVol  = numpy.log10(volume/volume.sum())
	lgDensSelect = lgDens[lgDens > 1.0]
	lgVolSelect = lgVol[lgDens > 1.0]

	lgVolSelect = lgVolSelect[lgDensSelect < 3.0]
	lgDensSelect = lgDensSelect[lgDensSelect < 3.0]
	fit= numpy.polyfit(lgDensSelect,lgVolSelect,1)
	scaling = fit[0]
	norm = 1e1**fit[1]
	print "scaling = " + str(scaling) + " norm = " + str( norm)
	print "fit = " + str(fit)
	dselect = 1e1**lgDensSelect 
	pl.plot(dselect, norm*(dselect)**scaling,'.',color=color,label="$\\alpha={0:3.2f}$".format(scaling))

pf = load("BB_hdf5_plt_cnt_0122")

# 512 file 122
posx = 7.49124383e+18
posy = 3.67059149e+19
posz = 3.23281834e+19

plotPDF(pf, posx, posy, posz, rmax = 1.0, label="1 pc",color='blue')
plotPDF(pf, posx, posy, posz, rmax = 2.0, label="2 pc",color='red')
plotPDF(pf, posx, posy, posz, rmax = 5.0, label="5 pc",color='black')

pl.xscale('log')
pl.yscale('log')
pl.ylim(1e-8, 1e-1)
pl.legend(loc='best')
pl.xlabel("$\\rho/\\rho_0$")
pl.ylabel("PDF")
pl.savefig("prob.pdf")
