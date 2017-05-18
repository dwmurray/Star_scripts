import matplotlib
matplotlib.use("Agg")
from yt.mods import *
import numpy
import math
import matplotlib.pyplot as p

rmax = 4. 
nbins = 100
pc = 3e18
meanDens = 3e-22
Msun = 2e33

pf = load("BB_hdf5_plt_cnt_0144")
dd = pf.h.all_data()
posxArr = dd["particle_posx"]
posyArr = dd["particle_posy"]
poszArr = dd["particle_posz"]
pmassArr = dd["particle_mass"]/Msun

posxArrSelect = posxArr
posyArrSelect = posyArr
poszArrSelect = poszArr
pmassArrSelect = pmassArr

avgDensBin = numpy.zeros(nbins)
rbin = (numpy.array(range(0,nbins))*1.0 + 1.0)*rmax/nbins
bins = 0 
for i in range(0,posxArrSelect.size) :
	posx = posxArrSelect[i]
	posy = posyArrSelect[i]
	posz = poszArrSelect[i]
	pmass = pmassArrSelect[i]

	sp = pf.h.sphere([posx,posy,posz], rmax/pf['pc'])
	
	cellVol = sp["CellVolume"]
	density = sp["Density"]
	x = sp["x"] - posx
	y = sp["y"] - posy
	z = sp["z"] - posz
	
	r = numpy.sqrt(x*x + y*y + z*z)


        volBin = numpy.zeros(nbins)
        densBin = numpy.zeros(nbins)
        for j in range(0,r.size) : 
		index = int( r[j]/(rmax*pc)*nbins)  
		if( index >= 0 and index < nbins) :
			volBin[index] += cellVol[j]
			densBin[index] += density[j]*cellVol[j]

	densBin = densBin/volBin
	maxDens = densBin[0]
	if(maxDens > 1000.*meanDens)  :
		densBin = densBin/maxDens
		avgDensBin = avgDensBin + densBin	
		bins=bins+1
#	p.plot(rbin, densBin, label="$m_* = {0:2.1f}$".format(pmass))

avgDensBin = avgDensBin/bins	
p.plot(rbin, avgDensBin)
numpy.savetxt( "test.out", zip(rbin,avgDensBin), fmt="%12.7e")
p.xlabel("r [pc]")
p.ylabel("$\\rho/\\rho_{\\rm max}$")

p.xlim(0.02, 5.0)
p.ylim(1e-5, 1.1)	

r = numpy.arange(0.2,0.7,0.1)
p.plot(r,0.7*(r/0.1)**-2.0,label="$r^{-2}$")
p.plot(r,1e-2*(r/0.1)**-3.0,label="$r^{-3}$")
p.xscale('log')
p.yscale('log')
p.legend(loc="best")
	
p.savefig("test.pdf")

