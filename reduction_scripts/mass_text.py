import matplotlib
matplotlib.use("Agg")
from matplotlib import rc
rc("text", usetex=True)
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
rcParams['axes.linewidth'] = 2
rcParams['lines.linewidth'] = 2
rcParams['xtick.major.width'] = 3    # major tick width in points
rcParams['xtick.major.size'] = 14    # major tick width in points                                         
rcParams['xtick.minor.size'] = 8    # major tick width in points                                                                        
rcParams['xtick.minor.width'] = 2    # major tick width in points                                                                       
rcParams['ytick.major.width'] = 3    # major tick width in points                                                                       
rcParams['ytick.minor.width'] = 2    # major tick width in points                                                                       
rcParams['ytick.major.size'] = 14    # major tick width in points
rcParams['ytick.minor.size'] = 8    # major tick width in points 

import h5py
import numpy as np
import os
import sys
import argparse
import glob
import matplotlib.pyplot as pl
import scipy.interpolate as si


def plotter(label="test", ls='-'):
    #pl.clf()
    pl.rc('text', usetex=True)
    pl.rc('font', family='serif')
    pl.loglog( t, M, label=label, ls=ls)#, 'b')#, label='$u_r$')

    pl.xticks(fontsize=25)
    pl.yticks(fontsize=25)

    pl.ylabel('$M_*$ ($M_\odot$)', fontsize = 25)
    pl.xlabel( "$t-t_*$ (${\\rm yrs}$)", fontsize = 25)
    pl.loglog( np.arange(2, 100, 0.1)*1e4, 1.0e-2 *(np.arange(2,100,0.1))**2, ls="--")
    #pl.ylim(1.0e-7, 1.0e0)
    #pl.xlim(1.0e4, 3.0e6)

    pl.rc('text', usetex=False)






all = np.loadtxt('fed_pad.txt')
t = []
M = []
for item in range(len(all)):
    t.append(all[item][0])
    M.append(all[item][1])
#    print t, M

plotter(label="checks", ls="-.")

all = np.loadtxt('pad.txt')
t = []
M = []
for item in range(len(all)):
    t.append(all[item][0])
    M.append(all[item][1])
#    print t, M


plotter(label="no checks")

pl.legend(loc='best')
pl.savefig('test_text.pdf')
