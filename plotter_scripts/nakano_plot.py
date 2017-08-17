import numpy as np
import matplotlib.pyplot as pl
from scipy.interpolate import UnivariateSpline

filein = 'ZAMS_ezer67.txt'

mass, radius, luminosity, age = np.loadtxt(filein, comments='#', unpack=True)
#print mass, radius, luminosity, age
#luminosity = mass**3.
#spl = UnivariateSpline(x, y)
spl = UnivariateSpline(mass, luminosity, k=4, s=0)

print mass, luminosity
lum_deriv = luminosity
#lum_deriv[0] = spl.derivatives(mass[0])[1]


xs = np.linspace(0.5, 100., 1000)
#xs = np.logspace(0.5, 100., 1000)

xs_deriv = np.zeros(len(mass))
#xs_deriv = np.zeros(len(xs))
print len(xs)
for i in range(len(mass)):
#    print i, xs[i]
    xs_deriv[i] = spl.derivatives(mass[i])[1]
#    xs_deriv[i] = spl.derivatives(xs[i])[1]
    print mass[i], luminosity[i], xs_deriv[i], xs_deriv[i] * mass[i] / luminosity[i]

#log_lum = np.zeros(len(mass))
#for i in range(len(mass)):#(0,1,2):#range(len(mass)):
#    lum_deriv[i] = spl.derivatives(mass[i])[1]
#    log_lum[i] = lum_deriv[i] / (mass[i]*mass[i])
#    print i, mass[i], luminosity[i], lum_deriv[i], log_lum[i]
#
#print lum_orig
#luminosity = lum_orig
#print mass, luminosity, lum_deriv, log_lum
#print len(xs), len(xs_deriv)

pl.clf()
pl.loglog(mass, luminosity, 'x', color='black')
#pl.loglog(mass, spl(mass), color='blue')
#pl.loglog(xs, spl(xs), color='blue')
pl.loglog(mass, xs_deriv, '^-', color='blue')
#pl.loglog(xs, xs_deriv, '^-', color='blue')
pl.xlim(3.e-1, 1.2e2)
pl.ylim(5.e-3, 2.e6)
pl.ylabel(r'Luminosity ($L_* / L_\odot$)')
pl.xlabel(r'Mass ($M_\odot$)')
pl.savefig('luminosity.pdf')

pl.clf()
pl.loglog(mass, radius, 'x-', color='black')
pl.xlim(3.e-1, 1.2e2)
pl.ylim(3.e-1, 2.e1)
pl.ylabel(r'Radius ($R_* / R_\odot$)')
pl.xlabel(r'Mass ($M_\odot$)')
pl.savefig('radius.pdf')
