import matplotlib
matplotlib.use('Agg')
from yt.mods import *
sys.path.append('./scripts/joishi/fourier_tools/fourier_tools/')

import matplotlib.pyplot as pl
import fourier_filter as ff
import numpy as np
import time
import math 

pf = load("BB_hdf5_plt_cnt_0080")
N = 512
cube = pf.h.covering_grid(level=5,left_edge=pf.domain_left_edge,right_edge=pf.domain_right_edge,dims=[N,N,N])
startime = time.time()
print 'start vel ' + str(time.time()-startime)
vx = cube["x-velocity"]
vy = cube["y-velocity"]
vz = cube["z-velocity"]
#magp = cube["magp"]
#magx = cube["magx"]
#magy = cube["magx"]
#magz = cube["magx"]
dens = cube["Density"]
#vel = cube["x-velocity"]


print 'start filter ' + str(time.time()-startime)
ffdat = ff.FourierFilter(vx)

print 'end filter ' + str(time.time()-startime)
bins = ffdat.bins[:ffdat.nx]

print 'start power ' + str(time.time()-startime)
fft_vx = na.fft.fftn(vx)/vx.size
fft_vy = na.fft.fftn(vy)/vy.size
fft_vz = na.fft.fftn(vz)/vz.size
print 'end power' + str(time.time()-startime)
divspec = na.zeros(ffdat.nx)
curlspec = na.zeros(ffdat.nx)
for i in range(ffdat.nx):
    for j in range(ffdat.nx) :
        for k in range(ffdat.nx) :
            bin = int(math.floor(math.sqrt(i*i + j*j + k*k)))
            if(bin < ffdat.nx) : 
                divspec[bin] = divspec[bin] + abs(i*fft_vx[i,j,k] + j*fft_vy[i,j,k] + k*fft_vy[i,j,k])**2
                curlspec[bin] = curlspec[bin] + abs(j*fft_vz[i,j,k] - k*fft_vy[i,j,k])**2 + abs(k*fft_vx[i,j,k] - i*fft_vy[i,j,k])**2 + abs(i*fft_vy[i,j,k] - j*fft_vx[i,j,k])**2

rho = 3e-22
L = 30e19
M = rho*L*L*L
Grav = 6.67e-8
cs = 0.265e5
#meanB = math.sqrt(8.*3.141*magp.mean())
#meanVa2 = (magp/dens).mean()
#print 'mean field = ' + str(meanB)
#meanrho = dens.mean()
#print 'mean density = ' + str(meanrho)
#print 'mean beta = ' + str(cs*cs/(meanVa2))
#print 'mean beta: 2 = ' + str(cs*cs/(meanB*meanB/(4.*3.141*meanrho)))
#print 'max b field = ' + str(math.sqrt(4.*3.141*(magx*magx + magy*magy + magz*magz).max() ))
#spec = na.array([power[ffdat.get_shell(bin)].sum() for bin in range(ffdat.nx)])
print 'end spec ' + str(time.time()-startime)

pl.plot(bins, divspec)
pl.plot(bins, curlspec, ls=':')
pl.xlim(1e-3,1.)
pl.yscale('log')
pl.xscale('log')
pl.xlabel('k')
pl.ylabel('$k^{2}P_{v}$')
pl.savefig('test.pdf')
numpy.savetxt("total_solenoid.data", zip(bins,divspec, curlspec))
