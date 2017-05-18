import matplotlib
matplotlib.use('Agg')
from yt.mods import *
sys.path.append('./scripts/joishi/fourier_tools/fourier_tools/')

import matplotlib.pyplot as pl
import fourier_filter as ff
import numpy as np
import time
import math 

pf = load("BB_hdf5_plt_cnt_0100")
N = 256 
cube = pf.h.covering_grid(level=4,left_edge=pf.domain_left_edge,right_edge=pf.domain_right_edge,dims=[N,N,N])
startime = time.time()
print 'start vel ' + str(time.time()-startime)
vel = cube["VelocityMagnitude"]
magp = cube["MagneticPressure"]
#magx = cube["magx"]
#magy = cube["magx"]
#magz = cube["magx"]
dens = cube["Density"]
#vel = cube["x-velocity"]


print 'start filter ' + str(time.time()-startime)
ffdat = ff.FourierFilter(vel)

print 'end filter ' + str(time.time()-startime)
bins = ffdat.bins[:ffdat.nx]

print 'start power ' + str(time.time()-startime)
power = na.abs(na.fft.fftn(vel)/vel.size)**2
print 'end power' + str(time.time()-startime)
spec = na.zeros(ffdat.nx)
for i in range(ffdat.nx):
    for j in range(ffdat.nx) :
        for k in range(ffdat.nx) :
            bin = int(math.floor(math.sqrt(i*i + j*j + k*k)))
            if(bin < ffdat.nx) : 
                spec[bin] = spec[bin] + power[i,j,k]

sigma = math.sqrt((vel*vel).mean())
rho = 3e-22
L = 30e19
M = rho*L*L*L
Grav = 6.67e-8
cs = 0.265e5
alpha_vir = 5.*sigma*sigma * L /(6.*Grav*M) # from Padoan and Norlund in the text pg. 7
print 'alpha_vir = ' + str(alpha_vir)
print 'mach number = ' + str(sigma/cs)
meanB = np.sqrt(8.*3.141*magp).mean()
#meanVa2 = (magp/dens).mean()
print 'mean field = ' + str(meanB)
#meanrho = dens.mean()
#print 'mean density = ' + str(meanrho)
#print 'mean beta = ' + str(cs*cs/(meanVa2))
#print 'mean beta: 2 = ' + str(cs*cs/(meanB*meanB/(4.*3.141*meanrho)))
print 'max b field = ' + str(math.sqrt(8.*3.141*magp.max()))
norm = (spec*bins*bins).sum()
#spec = na.array([power[ffdat.get_shell(bin)].sum() for bin in range(ffdat.nx)])
print 'end spec ' + str(time.time()-startime)

pl.plot(bins, spec*bins*bins/norm)
pl.xlim(1e-3,1.)
pl.yscale('log')
pl.xscale('log')
pl.xlabel('k')
pl.ylabel('$k^{2}P_{v}$')
pl.savefig('test.pdf')
numpy.savetxt("total_power_spec.data", zip(bins,spec))
