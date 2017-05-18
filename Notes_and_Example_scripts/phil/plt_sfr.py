import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as pl
import numpy

arr = numpy.loadtxt( "test.out")

index = arr[:,0]
time = arr[:,1]
t0   = arr[0,0]
star = arr[:,2]

print star

#time = 3.15e12*(index - 50)
t0 = time[0]

rho_mean = 3e-22
L = 4.8e19
Msun = 1.9885e33
total_mass = rho_mean*L*L*L/Msun
tff  = (3.*3.141/(32.*6.67e-8*rho_mean))**0.5

pl.plot((time-t0)/tff,star/total_mass)
pl.savefig("test.pdf")
