"""
Plot Vpow of 256 and 512
"""

import numpy as np
import scipy.optimize as so
import scipy.integrate as si
from matplotlib import rc
import matplotlib.pyplot as pl
import pdb

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

#Load data
hddir = '/mnt/raid-project/murray/elee/enzo/ic_pn/'
mhddir = '/mnt/raid-project/murray/elee/enzo/mhd_ic_pn_0.49microG/'

bins_256, spec_256 = np.loadtxt(hddir+'res_256/vpow_0000_0050_v.txt', unpack=True)
bins_512, spec_512 = np.loadtxt(hddir+'res_512/vpow_0000_0066_v.txt', unpack=True)
bins_m256, spec_m256 = np.loadtxt(mhddir+'res_256/HLL_PLM_ConsRecOff/vpow_0000_0135_v.txt', unpack=True)

#Normalize the power spectrum
norm_256 = si.simps(spec_256, x=bins_256)
norm_512 = si.simps(spec_512, x=bins_512)
norm_m256 = si.simps(spec_m256, x=bins_m256)

#Fit a power law
lfit = lambda p, x: p[0]*x + p[1]
residual = lambda p, x, f: lfit(p, x) - f

fiti_256 = np.where((bins_256 <= 1e-1)*(bins_256 > 0) > 0)[0]
fiti_512 = np.where((bins_512 <= 1e-1)*(bins_512 > 0) > 0)[0]
fiti_m256 = np.where((bins_m256 <= 1e-1)*(bins_m256 > 0) > 0)[0]

coeff_256 = so.leastsq(residual, np.array([-2., 0.]), args=(np.log10(bins_256)[fiti_256], np.log10(spec_256)[fiti_256]), full_output=True)
coeff_512 = so.leastsq(residual, np.array([-2., 0.]), args=(np.log10(bins_512)[fiti_512], np.log10(spec_512)[fiti_512]), full_output=True)
coeff_m256 = so.leastsq(residual, np.array([-2., 0.]), args=(np.log10(bins_m256)[fiti_m256], np.log10(spec_m256)[fiti_m256]), full_output=True)

success_256 = False
print 'HD256'
if coeff_256[-1] < 1 or coeff_256[-1] > 4:
    print 'FIT FAILED!'
    pdb.set_trace()
else:
    print 'FITTED!'
    success_256 = True
    print 'Power = %4.2f +/- %4.2f' %(coeff_256[0][0], coeff_256[1][0][0]**0.5)
print ''

success_512 = False
print 'HD512'
if coeff_512[-1] < 1 or coeff_512[-1] > 4:
    print 'FIT FAILED!'
    pdb.set_trace()
else:
    print 'FITTED!'
    success_512 = True
    print 'Power = %4.2f +/- %4.2f' %(coeff_512[0][0], coeff_512[1][0][0]**0.5)
print ''

success_m256 = False
print 'MHD256'
if coeff_m256[-1] < 1 or coeff_m256[-1] > 4:
    print 'FIT FAILED!'
    pdb.set_trace()
else:
    print 'FITTED!'
    success_m256 = True
    print 'Power = %4.2f +/- %4.2f' %(coeff_m256[0][0], coeff_m256[1][0][0]**0.5)
print ''

pl.clf()
pl.plot(bins_256, spec_256*bins_256**2./norm_256, 'k-')
pl.plot(bins_512, spec_512*bins_512**2./norm_512, 'b-')
pl.plot(bins_m256, spec_m256*bins_m256**2./norm_m256, 'r-')
if success_256: pl.plot(bins_256[1:], 10.**coeff_256[0][1]*bins_256[1:]**coeff_256[0][0]*bins_256[1:]**2./norm_256, 'k--')
if success_512: pl.plot(bins_512[1:], 10.**coeff_512[0][1]*bins_512[1:]**coeff_512[0][0]*bins_512[1:]**2./norm_512, 'b--')
if success_m256: pl.plot(bins_m256[1:], 10.**coeff_m256[0][1]*bins_m256[1:]**coeff_m256[0][0]*bins_m256[1:]**2./norm_m256, 'r--')
pl.yscale('log')
pl.xscale('log')
pl.xlabel('k')
#pl.ylim([1e-4, 3e-3])
#pl.xlim([3e-3, 1])
pl.ylabel(r'$k^{2}P_{v}$')
pl.show()






