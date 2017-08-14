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
from matplotlib import rcParams

parser = argparse.ArgumentParser(description = "start number to end number")

parser.add_argument('start', metavar='N1', type=int)
parser.add_argument('end', metavar='N2', type=int)
parser.add_argument('step', metavar='N3', type=int)
parser.add_argument('--one', action='store_true')
parser.add_argument('--two', action='store_true')
parser.add_argument('--plot', action='store_true')
args = parser.parse_args()


dirs = ["ultrares/jet", "ultrares/nojet", "highres/jet", "highres/nojet"]
labels = ["$32K^3$ jet", "$32K^3$ no jet", "$16K^3$ jet", "$16K^3$ no jet"]
ltypes = ["solid", "solid", "dashed", "dashed"]
lweights = [4, 2, 4, 2]
lcolors = ["b", "r", "b", "r"]

if args.one:
    dirs = ["ultrares/jet", "ultrares/nojet"]
    labels = ["$32K^3$ jet", "$32K^3$ no jet"]
    ltypes = ["solid", "solid"]
    lweights = [4, 2]
    lcolors = ["b", "r"]

if args.two:
    dirs = ["highres/jet", "highres/nojet"]
    labels = ["$16K^3$ jet", "$16K^3$ no jet"]
    ltypes = ["dashed", "dashed"]
    lweights = [4, 2]
    lcolors = ["b", "r"]


if not args.plot:
    print 'Got here.'
    for dir, label, ltype, lw, lc in zip(dirs,labels,ltypes,lweights,lcolors): 
        first_time = True
        t_start = 0.
        tarray = []
        vtot = []
        prefix = "output_"
        savefileout = "kin_E_{1}.txt".format(dir, label)
        print dir
        for i in range(args.start,args.end,args.step) :	
            infofile = "{0}/{1}".format(dir,prefix)+"{0:05d}/info_{0:05d}.txt".format(i)
            sinkfile = "{0}/{1}".format(dir,prefix)+"{0:05d}/sink_{0:05d}.info".format(i)
            if( not os.path.isfile(infofile)) : 
                continue
            if( not os.path.isfile(sinkfile)) : 
                continue
            print "opening {0}".format(infofile)
            v_tot = 0.
            try : 
                pf = yt.load(infofile)
                mass, xstar, ystar, zstar = np.loadtxt( sinkfile, usecols=[1,2,3,4], unpack=True, skiprows=3, comments="=")
            except ValueError : 
                print infofile + " empty or " + sinkfile + " empty" 
                continue

            print 'loading data'
            ad = pf.all_data()
            time = pf.current_time # Yt does not convert to seconds from code time automatically. # on home computer.
            time = time.in_cgs()
            time = time / 3.1558e7 # / 1e6 # Convert to years.# convert to millions of years

            cellMass = ad["cell_mass"]
            vx = ad["x-velocity"].in_cgs()
            vy = ad["y-velocity"].in_cgs()
            vz = ad["z-velocity"].in_cgs()

            cellMass = np.array(cellMass)
            vx = np.array(vx)
            vy = np.array(vy)
            vz = np.array(vz)
            
            total_mass = cellMass.sum()
            
            vx = vx / 1e5
            vy = vy / 1e5
            vz = vz / 1e5
            px = cellMass * vx
            py = cellMass * vy
            pz = cellMass * vz
            v_cmx = px / total_mass
            v_cmy = py / total_mass
            v_cmz = pz / total_mass
            
            vx = vx - v_cmx
            vy = vy - v_cmy
            vz = vz - v_cmz
            
            v_tot_sq = vx*vx + vy*vy + vz*vz
            
            p = cellMass * v_tot_sq
            total = p.sum()
            v_tot = total / total_mass
        
            if( first_time) :
                t_start = time
                print 'start'
                first_time = False
            if( not first_time) :
                print 'T'
                vtot.append(v_tot)
                tarray.append(time-t_start)
            print 'plotting tarray and vtot'
            print tarray, vtot
            print 'for dir, label:', dir, label
        np.savetxt(savefileout, zip(tarray, vtot), fmt="%15.9E")
        plt.plot(tarray, vtot, label="{0}".format(label), linewidth=lw, color=lc, ls=ltype)
        plt.savefig('{0}_plotted.pdf'.format(dir))

if args.plot:
    
    rcParams['axes.linewidth']    = 2
    rcParams['lines.linewidth']   = 2
    rcParams['xtick.major.width'] = 3    # major tick width in points
    rcParams['xtick.major.size']  = 14    # major tick width in points
    rcParams['xtick.minor.size']  = 8    # major tick width in points
    rcParams['xtick.minor.width'] = 2    # major tick width in points
    rcParams['ytick.major.width'] = 3    # major tick width in points
    rcParams['ytick.minor.width'] = 2    # major tick width in points
    rcParams['ytick.major.size']  = 14    # major tick width in points
    rcParams['ytick.minor.size']  = 8    # major tick width in points
    rcParams['xtick.labelsize']   = 25
    rcParams['ytick.labelsize']   = 25

    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    plt.minorticks_on()
    plt.tick_params('both',length=8, width=1, which='minor')
    plt.tick_params('both',length=10, width=1.5, which='major')
    plt.gcf().subplots_adjust(bottom=0.15)
    #plt.tight_layout()

    files = ["kin_E_32k_jet.txt", "kin_E_32k_nojet.txt", "kin_E_16k_jet.txt", "kin_E_16k_nojet.txt"]
    labels = ["$32K^3$ jet", "$32K^3$ no jet", "$16K^3$ jet", "$16K^3$ no jet"]
    ltypes = ["solid", "solid", "dashed", "dashed"]
    lweights = [4, 2, 4, 2]
    lcolors = ["b", "r", "b", "r"]
    for file, label, ltype, lw, lc in zip (files, labels, ltypes, lweights, lcolors):
        tarray, vtot = np.loadtxt(file, unpack=True)
        plt.plot(tarray / 1.e6, vtot, label="{0}".format(label), linewidth=lw, color=lc, ls=ltype)
    plt.ylim(1e0, 1.2e1)
    plt.xlim(0e0, 1.1e0)
    plt.legend(loc='best')
    plt.ylabel(r'$<v_{rms}^2>$ $({\rm km^2\, s}^{-2})$', fontsize = 25)
    plt.xlabel(r'$t-t_{*}$ $({\rm Myrs})$', fontsize=25)
    plt.savefig('vsquared.pdf')
