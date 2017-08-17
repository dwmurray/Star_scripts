import numpy as np
import matplotlib.pyplot as pl
import math

G_N = 6.67e-8
Msun = 1.998e33
M_bh = 10.*Msun
c = 3e10
cs = 0.2*c
rs = 2.0*G_N*M_bh/(c*c)

# inner boundary conditions
rho0 = 0.1
vr0 = cs*1.0

# courant condition
courant_condition = 0.05

def RHS( r, rho, vr, P) :
    # r, rho, vr, P are all arrays of some size

    r2 = (r*r)
    N = r.size
    dr = r[1:N] - r[0:N-1]

    drhodt = np.zeros(N)
    dvrdt = np.zeros(N)

    drhodt[1:N] = -1.0/r2[1:N] * 1./dr*((rho*vr*r2)[1:N] - (rho*vr*r2)[0:N-1])
    dvrdt[1:N] = -vr[1:N]/dr*(vr[1:N] - vr[0:N-1]) - 1.0/rho[1:N]*1.0/dr*(P[1:N] - P[0:N-1]) - G_N*M_bh/r2[1:N]

    return [drhodt, dvrdt]

def applyBC( rho1, vr1) : 
    rho1[0] = rho0
    vr1[0] = vr0
    
    return [rho1, vr1]

def evolve( r, rho, vr, P, dt) : 
    
    # first evolve
    drhodt, dvrdt = RHS( r, rho, vr, P) 
    rho1 = rho + drhodt * dt
    vr1 = vr + dvrdt * dt
     
    # apply boundary conditions
    rho1, vr1 = applyBC( rho1, vr1)
 
    P1 = rho1*cs*cs   

    return rho1, vr1, P1

def cal_courant( r, vr) :
    dr = r[1:N] - r[0:N-1]
    dt = np.abs(dr/vr[1:N] * courant_condition)

    return np.min(dt)

# main program

# set up the logarithmic grid in r from 10 Rs to 10000 Rs
N = 400 # number of grid points
# inner and outer grid point
lgrin = math.log10( 10.*rs)
lgrout = math.log10( 1e4*rs) 
lgr = np.arange( lgrin, lgrout, (lgrout-lgrin)/N)

r = 1e1**lgr
rho = rho0/(r/r[0])**1.0 # for no gravity
rho = rho0*np.exp( G_N*M_bh/(r*cs*cs) - G_N*M_bh/(r[0]*cs*cs)) # equilibrium for gravity
vr = np.ones( N) * cs * 10.0

# locks in inner boundary condition
rho, vr = applyBC( rho, vr)

P = rho*cs*cs

t = 0 
tend = r[-1]/cs * 10. # 10x the sound crossing time
dt = r[0]/c*0.00001

i = 0
dtframe = r[0]/cs*0.1
while (t < tend) : 
    if( t > dtframe* i) : 
        pl.clf()
        pl.loglog( r/rs, vr/cs)
        pl.xlim(5.0,50.)
        pl.savefig( "frame_{0:04d}.png".format(i))
        i = i + 1    
    rho, vr, P = evolve( r, rho, vr, P, dt)  
    t = t + dt
    # figure out the next time step
    dt = cal_courant( r, vr)
    
    print t/tend, dtframe/tend, cs/(G_N*M_bh/r[0])**0.5
    






