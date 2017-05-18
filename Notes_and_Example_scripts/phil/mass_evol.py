from yt.mods import *
from scipy.optimize import leastsq
import pylab as pl
import sys
import pdb

# Accept command-line arguments
datadir = sys.argv[1]
pltdir = sys.argv[2]
start = int(sys.argv[3])
end = int(sys.argv[4])
step = int(sys.argv[5])

def models(t, m):
    """
    Return simple and fitted models to mass fraction evolution
    t: time array
    m: mass fraction array
    """

    t1 = (t-t[0])[2]
    t2 = (t-t[0])[len(t)-2]
    mu1 = m[2]
    mu2 = m[len(m)-2]
    
    N0 = ((t2*mu1**(-0.25) - t1*mu2**(-0.25))/(t2-t1))**(-4.)
    alpha = (mu2**(-0.25) - mu1**(-0.25))/(mu1**(-0.25)*t2 - mu2**(-0.25)*t1)
    simple_model = N0*(1+alpha*(t-t[0]))**(-4.)

    fp = lambda p, t_in: p[0]*(1+p[1]*t_in)**(-4.)
    err = lambda p, mu, t_in: fp(p, t_in)-mu

    p0 = [N0, alpha]
    p, success = leastsq(err, p0, args=(m, t-t[0]), maxfev=10000)

    if success: print "SUCCESS!"

    model = p[0]*(1+p[1]*(t-t[0]))**(-4.)

    return simple_model, model

def plot_frac(t, m, ind, t_dyn, log_yr, yr):
    """
    Plot maximum mass fraction evolution
    t: time array
    m: mass fraction array
    ind: indicating which maximum (1 for max, 2 for 2nd maximum, 3 for 3rd maximum)
    t_dyn: dynamical time
    log_yr: yrange in log-log plot
    yr: yrange in non-log plot
    """

    simple_model, model = models(t, m)

    #pdb.set_trace()

    pl.clf()
    pl.plot(na.log10((t-t[0])/t_dyn), na.log10(m))
    pl.plot(na.log10((t-t[0])/t_dyn), na.log10(simple_model), 'r--')
    pl.plot(na.log10((t-t[0])/t_dyn), na.log10(model), 'g--')
    pl.xlabel('time(Tdyn = %3.1f)' %t_dyn)
    pl.ylabel('%i Max Star Mass / Cloud Mass' %(ind))
    pl.ylim(log_yr)
    pl.savefig(pltdir+'max_%i_mass_evol_log.eps' %(ind))

    pl.clf()
    pl.plot((t-t[0])/t_dyn, m)
    pl.plot((t-t[0])/t_dyn, simple_model, 'r--')
    pl.plot((t-t[0])/t_dyn, model, 'g--')
    pl.xlabel('time(Tdyn = %3.1f)' %t_dyn)
    if ind < 4:
        pl.ylabel('%i Max Star Mass / Cloud Mass' %(ind))
    elif ind == 4:
        pl.ylabel('Cumul Max Star Mass / Cloud Mass')
    pl.ylim(yr)
    pl.savefig(pltdir+'max_%i_mass_evol.eps' %(ind))

numstararr = na.array(range(start, end, step))*0.#na.zeros((end-start)/step)
massarr = na.array(range(start, end, step))*0.#na.zeros((end-start)/step)
maxmarr = na.array(range(start, end, step))*0.#na.zeros((end-start)/step)
max2marr = na.array(range(start, end, step))*0.#na.zeros((end-start)/step)
max3marr = na.array(range(start, end, step))*0.#na.zeros((end-start)/step)
max4marr = na.array(range(start, end, step))*0.
max5marr = na.array(range(start, end, step))*0.
max6marr = na.array(range(start, end, step))*0.
max7marr = na.array(range(start, end, step))*0.

#tarr = na.array(range(start, end, step))*0.5#na.linspace(start, end-1, (end-start)/step)*0.5 #in Myr
tarr = na.array(range(start, end, step))*0.5
#tmp_tarr = na.copy(tarr)
msun = 2.e33
chk_mass_arr = na.array(range(start, end, step))*0.

#find the maximum mass star
pf = load(datadir+"DD%04i/data%04i" %(end-1, end-1))
end_particles = pf.h.find_particles_by_type(4)
end_particle_ind = end_particles["particle_index"][na.argsort(end_particles["particle_mass"])]

maxi = -1
max2i = -1
max3i = -1
max4i = -1
max5i = -1
max6i = -1
max7i = -1

if len(end_particle_ind) > 0:
    maxi = end_particle_ind[len(end_particle_ind)-1]
if len(end_particle_ind) > 1:
    max2i = end_particle_ind[len(end_particle_ind)-2]
if len(end_particle_ind) > 2:
    max3i = end_particle_ind[len(end_particle_ind)-3]
if len(end_particle_ind) > 3:
    max4i = end_particle_ind[len(end_particle_ind)-4]
if len(end_particle_ind) > 4:
    max5i = end_particle_ind[len(end_particle_ind)-5]
if len(end_particle_ind) > 5:
    max6i = end_particle_ind[len(end_particle_ind)-6]
if len(end_particle_ind) > 6:
    max7i = end_particle_ind[len(end_particle_ind)-7]

for i in range(start, end, step):
    pf = load(datadir+"DD%04i/data%04i" %(i, i))
    sp = pf.h.sphere([0.5, 0.5, 0.5], 1.0)
    
    if i == start:
        sp = pf.h.sphere([0.5, 0.5, 0.5], 1.0)
        total_mass = sp["CellMassMsun"].sum()
  	t_dyn = pf.get_parameter('LengthUnits')/sp["VelocityMagnitude"].mean()/3600./24./365.25/1e6
    	#t_dyn = 10.

    tarr[(i-start)/step] = pf["InitialTime"]#/t_dyn
    particles = pf.h.find_particles_by_type(4)
    chk_mass_arr[(i-start)/step] = sp["CellMassMsun"].sum()

    if type(particles) != types.NoneType:
        chk_mass_arr[(i-start)/step] += sp["ParticleMassMsun"].sum()
        index = particles["particle_index"]
        mass = particles['particle_mass']*\
            ((1./pf["TopGridDimensions"][0])**3)*pf.get_parameter('DensityUnits')*(pf.get_parameter('LengthUnits')**3)
        numstararr[(i-start)/step] = len(mass)
        massarr[(i-start)/step] = mass.sum()

        #pdb.set_trace()

        if index.max() >= maxi and maxi > -1:
            maxmarr[(i-start)/step] = mass[na.where(index == maxi)].max()
        if index.max() >= max2i and max2i > -1:
            max2marr[(i-start)/step] = mass[na.where(index == max2i)].max()
        if index.max() >= max3i and max3i > -1:
            max3marr[(i-start)/step] = mass[na.where(index == max3i)].max()
        if index.max() >= max4i and max4i > -1:
            max4marr[(i-start)/step] = mass[na.where(index == max4i)].max()
        if index.max() >= max5i and max5i > -1:
            max5marr[(i-start)/step] = mass[na.where(index == max5i)].max()
        if index.max() >= max6i and max6i > -1:
            max6marr[(i-start)/step] = mass[na.where(index == max6i)].max()
        if index.max() >= max7i and max7i > -1:
            max7marr[(i-start)/step] = mass[na.where(index == max7i)].max()

        massarr[(i-start)/step] = (massarr[(i-start)/step]/msun)/total_mass
        maxmarr[(i-start)/step] = (maxmarr[(i-start)/step]/msun)/total_mass
        max2marr[(i-start)/step] = (max2marr[(i-start)/step]/msun)/total_mass
        max3marr[(i-start)/step] = (max3marr[(i-start)/step]/msun)/total_mass
        max4marr[(i-start)/step] = (max4marr[(i-start)/step]/msun)/total_mass
        max5marr[(i-start)/step] = (max5marr[(i-start)/step]/msun)/total_mass
        max6marr[(i-start)/step] = (max6marr[(i-start)/step]/msun)/total_mass
        max7marr[(i-start)/step] = (max7marr[(i-start)/step]/msun)/total_mass

        pl.clf()
        pl.hist(mass/msun,bins=20,log=True)
        pl.xlabel('Mass (Msun)')
        pl.ylabel('Number')
        pl.savefig(pltdir+'part_mass_hist_%04i.eps' %(i))
        
pl.plot(tarr, na.log10(chk_mass_arr))
pl.plot(tarr, na.log10(na.ones(len(tarr))*total_mass))
pl.xlabel('time(Myr)')
pl.ylabel('Log(Total Mass in Msun)')
pl.ylim([4.26, 4.28])
pl.savefig(pltdir+'chk_total_mass.eps')

pl.clf()
pl.plot((tarr-tarr[0])/t_dyn, massarr)
pl.xlabel('time(Tdyn = %3.1f)' %t_dyn)
pl.ylabel('M*/M')
pl.savefig(pltdir+'star_massf_evol.eps')

na.savetxt(pltdir+'star_massf_evol.txt', zip((tarr-tarr[0])/t_dyn, massarr))

pl.clf()
pl.plot(tarr, numstararr)
pl.xlabel('time(Myr)')
pl.ylabel('Number of Stars')
pl.savefig(pltdir+'n_star_evol.eps')

"""
#t_dyn = 4.8
log_yr = na.log10([maxmarr.min(), maxmarr.max()])+[-0.01, 0.01] #[-3.0, 0.0]
#log_yr = [0.0, 0.5] 
yr = [0.0, maxmarr.max()*1.5]
#yr = [0.0, 0.5]

plot_frac(tarr, maxmarr, 1, t_dyn, log_yr, yr)
plot_frac(tarr, max2marr, 2, t_dyn, log_yr, yr)
plot_frac(tarr, max3marr, 3, t_dyn, log_yr, yr)
plot_frac(tarr, maxmarr+max2marr+max3marr, 4, t_dyn, log_yr, yr)
"""

na.savetxt(pltdir+'indiv_star_massf_evol.txt', zip((tarr-tarr[0])/t_dyn, maxmarr, max2marr, max3marr, max4marr, max5marr, max6marr, max7marr))
