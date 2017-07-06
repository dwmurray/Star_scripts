import shutil
import matplotlib 
matplotlib.use("Agg")
import argparse
import os
import sys
import yt
import numpy as np
import matplotlib.pyplot as plt
import math

parser = argparse.ArgumentParser(description = "start number to end number")

parser.add_argument('start', metavar='N1', type=int)
parser.add_argument('end', metavar='N2', type=int)
parser.add_argument('step', metavar='N3', type=int)
parser.add_argument('--kink', action='store_true')
args = parser.parse_args()


dirs = ["ultrares/jet", "ultrares/nojet", "highres/jet", "highres/nojet"]
labels = ["$32K^3$ jet", "no jet", "$16K^3$ jet", "no jet"]
ltypes = ["solid", "dashed", "solid", "dashed"]
lweights = [4, 4, 2, 2]



for dir, label, ltype, lw in zip(dirs,labels,ltypes,lweights): 
    first_time = True
    t_start = 0.
    tarray = []
    vtot = []
    prefix = "output_"
    for i in range(args.start,args.end,args.step) :	
        file = "{0}/{1}".format(dir,prefix)+"{0:05d}/info_{0:05d}.txt".format(i)
        sinkfile = "{0}/{1}".format(dir,prefix)+"{0:05d}/sink_{0:05d}.info".format(i)
        if( not os.path.isfile(file)) : 
            continue
        if( not os.path.isfile(sinkfile)) : 
            continue
        print "opening {0}".format(sinkfile)
        vtot = None
        try : 
            #Particle_ID_list, partMass, r_star, xstar, ystar, zstar, vxstar, vystar, vzstar, poly_n, md, polystate, pjet = np.loadtxt( sinkfile, unpack=True, skiprows=3, comments="=")
                    #mass, xstar, ystar, zstar = numpy.loadtxt( sinkfile, usecols=[1,2,3,4], unpack=True, skiprows=3, comments="=")
            pf = yt.load(file)
            ad = pf.all_data()
            time = pf.current_time # Yt converts to seconds from code time automatically.
            cellMass = ad["cell_mass"]
            total_mass = cellMass.sum()
            
            vx = ad["x-velocity"].in_cgs()
            vy = ad["y-velocity"].in_cgs()
            vz = ad["z-velocity"].in_cgs()
            v_tot_sq = vx*vx + vy*vy + vz*vz

            p = cellMass * v_tot_sq
            total = p.sum()
            vtot = total / total_mass

            if( first_time) :
                t_start = time
                first_time = False
            if( not first_time) :
                vtot.append(v_tot)
                tarray.append(time-t_start)
        except ValueError : 
            print infofile + " empty" 

#        gasMass = 4.8e19**3 * 3e-22 / 2e33
        tarray = tarray / 3.1558e7 # Convert to years.
        tarray = tarray / 1e6 # convert to millions of years
	#tarray = numpy.log10(tarray) 
	#tarray = 1e1**tarray/3.1558e7
	#mtot = numpy.log10(mtot) - math.log10(gasMass)
	#mtot = 1e1**mtot

	for t, v in zip(tarray, vtot) :
		print t, v

        plt.plot(tarray, vtot, label="{0}".format(label))

plt.ylim(1e0, 1.2e1)
plt.xlim(0e0, 3e0)
plt.ylabel(r'$<v_{rms}>$ $({\rm km\, s}^{-1})$', fontsize = 25)
plt.xlabel(r'$t$ $({\rm Myrs})$', fontsize=25)
	
#	lv = numpy.log10(vtot)
#        lt = numpy.log10(tarray)
#
#	lt_min = 5.5
#	lt_max = 6.1
#        p = numpy.polyfit( lt[lt>lt_min], lm[lt>lt_min], 1)


#
#cellMass = ad["cell_mass"]
#total_mass = cellMass.sum()
#
#vx = ad["x-velocity"].in_cgs()
#vy = ad["y-velocity"].in_cgs()
#vz = ad["z-velocity"].in_cgs()
#v_tot_sq = vx*vx + vy*vy + vz*vz
#
#p = cellMass * v_tot_sq
#total = p.sum()
#
#v_tot = total / total_mass
#v_tot / 1e10
#
