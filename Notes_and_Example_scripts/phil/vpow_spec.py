"""
Program to plot velocity power spectrum
at the given Mach number
--designed to check expected scaling law
at the fully-turbulent regime

@author: Eve Lee
Jul 5th 2012
"""

from yt.mods import *
import scipy.optimize as so
import scipy.integrate as si
from matplotlib import rc
import matplotlib.pyplot as pl
import pow_spec
import pdb
import sys
sys.path.append('./joishi/fourier_tools/fourier_tools/')
import fourier_filter as ff

## Set up plot parameters
rc('font', size=20, **{'family':'serif'})
rc('text', usetex=True)
rc('xtick', labelsize=15)
rc('xtick.major', size=8)
rc('xtick.minor', size=4)
rc('ytick', labelsize=15)
rc('ytick.major', size=8)
rc('ytick.minor', size=4)
rc('axes', labelsize=20)
rc('legend', fontsize=10)

#Read in system argument
datadir = sys.argv[1]
outdir = sys.argv[2]
start = sys.argv[3]
end = sys.argv[4]
Mach = float(sys.argv[5])

#Load data (should be time, Mach number, mass-weighted Mach number)
t, mach, wmach = na.loadtxt(datadir+'mach_'+start+'_'+end+'.txt', unpack=True)
t, vir, vir_class, vir_enzo = na.loadtxt(datadir+'vir_'+start+'_'+end+'.txt', unpack=True)

#Compare mach and wmach evolution
pl.clf()
pl.plot(t, mach)
pl.plot(t, wmach)
pl.legend(('Mach', 'wMach'), loc='best')
pl.xlabel('Time (Myr)')
pl.show()

#Compare different virial parameter
pl.clf()
pl.plot(t, vir)
pl.plot(t, vir_class)
pl.plot(t, vir_enzo)
pl.legend(('KE/PE','Classical','Enzo'), loc='best')
pl.ylabel(r'$\alpha_{vir}$')
pl.xlabel('Time (Myr)')
pl.show()

#Find the supposedly fully-turbulent data cube
goodi = na.argmin(na.abs(wmach - Mach))
pf = load(datadir+'DD%04i/data%04i' %(goodi, goodi))
grid = pf.h.grids[0]
res = grid.ActiveDimensions[0]
grid = 0. # to save memory
dim = pf.domain_dimensions
cube = pf.h.covering_grid(0, left_edge=[0.,0.,0.], right_edge=[1.,1.,1.], dims=dim, fields=["VelocityMagnitude"])
vel = cube["VelocityMagnitude"]
vx = cube["x-velocity"] - na.mean(cube["x-velocity"])
vy = cube["y-velocity"] - na.mean(cube["y-velocity"])
vz = cube["z-velocity"] - na.mean(cube["z-velocity"])
rho = cube["Density"]

#Calculate power spectrum
#kpplot, pplot = pow_spec.powspec(vel, 50)
ffdat = ff.FourierFilter(vel)
power = na.abs(na.fft.fftn(vel)/vel.size)**2.
spec = na.array([power[ffdat.get_shell(bin)].sum() for bin in range(ffdat.nx)])
bins = ffdat.bins[:ffdat.nx]
norm = si.simps(spec, x=bins)

ffdat_x = ff.FourierFilter(vx)
ffdat_y = ff.FourierFilter(vy)
ffdat_z = ff.FourierFilter(vz)

power2 = na.abs(na.fft.fftn(vx)/vx.size)**2. + \
    na.abs(na.fft.fftn(vy)/vy.size)**2. + \
    na.abs(na.fft.fftn(vz)/vz.size)**2.
spec2 = na.array([power2[ffdat_x.get_shell(bin)].sum() for bin in range(ffdat_x.nx)])
bins2 = ffdat_x.bins[:ffdat_x.nx]
norm2 = si.simps(spec2, x=bins2)

ffdat_rx = ff.FourierFilter(vx*(rho**(1./3.)))
ffdat_ry = ff.FourierFilter(vy*(rho**(1./3.)))
ffdat_rz = ff.FourierFilter(vz*(rho**(1./3.)))

power_r = na.abs(na.fft.fftn(vx*(rho**(1./3.)))/vx.size)**2. + \
    na.abs(na.fft.fftn(vy*(rho**(1./3.)))/vy.size)**2. + \
    na.abs(na.fft.fftn(vz*(rho**(1./3.)))/vz.size)**2.
spec_r = na.array([power_r[ffdat_rx.get_shell(bin)].sum() for bin in range(ffdat_rx.nx)])
bins_r = ffdat_rx.bins[:ffdat_rx.nx]
normr = si.simps(spec_r, x=bins_r)

print 'Mach = %4.2f Alpha_vir = %4.2f / %4.2f at Data No. %4i' %(wmach[goodi], vir[goodi], vir_class[goodi], goodi)

#Fit a power law
lfit = lambda p, x: p[0]*x + p[1]
residual = lambda p, x, f: lfit(p, x) - f

fiti = na.where((bins >= 1.e-2)*(bins <= 1e-1) > 0)[0]
fiti2 = na.where((bins2 >= 1.e-2)*(bins2 <= 1e-1) > 0)[0]
fiti_r = na.where((bins_r >= 1.e-2)*(bins_r <= 1e-1) > 0)[0]
#coeff = so.leastsq(residual, na.array([-4., 0.]), args=(na.log10(kpplot[fiti]), na.log10(pplot[fiti])), full_output=True)
coeff = so.leastsq(residual, na.array([-4., 0.]), args=(na.log10(bins)[fiti], na.log10(spec)[fiti]), full_output=True)
coeff2 = so.leastsq(residual, na.array([-4., 0.]), args=(na.log10(bins2)[fiti2], na.log10(spec2)[fiti2]), full_output=True)
coeff_r = so.leastsq(residual, na.array([-4., 0.]), args=(na.log10(bins_r)[fiti_r], na.log10(spec_r)[fiti_r]), full_output=True)

#Report fit result
success = False
print 'FFT(sigma)'
if coeff[-1] < 1 or coeff[-1] > 4:
    print 'FIT FAILED!'
    pdb.set_trace()
else:
    print 'FITTED!'
    success = True
    print 'Power = %4.2f +/- %4.2f' %(coeff[0][0], coeff[1][0][0]**0.5)
print ''

success2 = False
print 'FFT(v-mean_v)'
if coeff2[-1] < 1 or coeff2[-1] > 4:
    print 'FIT FAILED!'
    pdb.set_trace()
else:
    print 'FITTED!'
    success2 = True
    print 'Power = %4.2f +/- %4.2f' %(coeff2[0][0], coeff2[1][0][0]**0.5)
print ''

success_r = False
print 'FFT(rho-weighted v-mean_v)'
if coeff_r[-1] < 1 or coeff_r[-1] > 4:
    print 'FIT FAILED!'
    pdb.set_trace()
else:
    print 'FITTED!'
    success_r = True
    print 'Power = %4.2f +/- %4.2f' %(coeff_r[0][0], coeff_r[1][0][0]**0.5)
print ''

na.savetxt(outdir+'vpow_%04i_%04i_sigma.txt' %(int(start), int(end)), zip(bins, spec))
na.savetxt(outdir+'vpow_%04i_%04i_v.txt' %(int(start), int(end)), zip(bins2, spec2))
na.savetxt(outdir+'vpow_%04i_%04i_rhowv.txt' %(int(start), int(end)), zip(bins_r, spec_r))

pl.clf()
#pl.plot(kpplot, pplot)
pl.plot(bins, spec*bins**2./norm)
pl.plot(bins2, spec2*bins2**2./norm2)
pl.plot(bins_r, spec_r*bins_r**2./normr) #originally spec_r*1e6
#if success: pl.plot(kpplot, 10.**coeff[0][1]*kpplot**coeff[0][0], '--')

if success: pl.plot(bins, 10.**coeff[0][1]*bins**coeff[0][0]*bins**2./norm, '--')
if success2: pl.plot(bins2, 10.**coeff2[0][1]*bins2**coeff2[0][0]*bins2**2./norm2, '--')
if success_r: pl.plot(bins_r, 10.**coeff_r[0][1]*bins_r**coeff_r[0][0]*bins_r**2./normr, '--')

pl.yscale('log')
pl.xscale('log')
pl.xlabel('k')
#pl.ylim([1e-3, 1e3])
#pl.xlim([1e-3, 1])
if success: pl.legend(('FFT(SIGMA)', 'FFT(V)', 'FFT(u)','Fit-sig','Fit-v','Fit-u'), loc='best')
pl.ylabel(r'$k^{2}P_{v}$')
pl.show()

pl.clf()
pl.plot(bins2, spec2*bins2**2./norm2, 'k-')
if success: pl.plot(bins2, 10.**coeff2[0][1]*bins2**coeff2[0][0]*bins2**2./norm2, 'k--')
pl.yscale('log')
pl.xscale('log')
pl.xlabel('k')
pl.ylim([1e-4, 3e-3])
pl.xlim([2e-3, 1])
pl.ylabel(r'$k^{2}P_{v}$')
pl.show()



