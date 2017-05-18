import matplotlib
matplotlib.use('Agg')
from yt.mods import *
sys.path.append('./scripts/joishi/fourier_tools/fourier_tools/')

import matplotlib.pyplot as pl
import fourier_filter as ff
import numpy as np
import time
import math 
import argparse

parser = argparse.ArgumentParser(description = "filename")

parser.add_argument('filename', metavar='N1')
args = parser.parse_args()

pf = load(args.filename)
data = pf.h.all_data()
vel = data["VelocityMagnitude"]
magp = data["magp"]
#magx = cube["magx"]
#magy = cube["magx"]
#magz = cube["magx"]
dens = data["Density"]
#vel = cube["x-velocity"]

sigma = math.sqrt((vel*vel).mean())
rho = 3e-22
L = 150e18
M = rho*L*L*L
Grav = 6.67e-8
cs = 0.265e5
alpha_vir = 5.*sigma*sigma * L /(6.*Grav*M) # from Padoan and Norlund in the text pg. 7
print 'alpha_vir = ' + str(alpha_vir)
print 'mach number = ' + str(sigma/cs)
meanB = math.sqrt(8.*3.141*magp.mean())
meanVa2 = (2.*magp/dens).mean()
print 'mean field = ' + str(meanB)
meanrho = dens.mean()
print 'mean density = ' + str(meanrho)
print 'mean beta = ' + str(rho*cs*cs/magp.mean())
print 'mean beta: 2 = ' + str(cs*cs/(meanB*meanB/(8.*3.141*meanrho)))
print 'max b field = ' + str(math.sqrt(8.*3.141*magp.max() ))
#spec = na.array([power[ffdat.get_shell(bin)].sum() for bin in range(ffdat.nx)])
