
import numpy as np
import scipy.integrate as si
import matplotlib.pyplot as pl 


mp = 1.6e-24
pi = np.pi
parsec = 3e18
sec_in_yr = np.pi* 1e7
Msun = 2e33
G = 6.67e-8
c_s = 2.65e4
cs = c_s
r_sun = 6.95949e10
sigma_b = 5.67e-5

pl.clf()
pl.rc('text', usetex=False)
pl.rc('font', family='serif')

i = ['nakano_10x_31', 'nakano_31', 'nakano_guess', 'offner']#[0.5, 1, 2, 5, 10, 20, 50, 100]
for j in range(len(i)):
#    Input = 'mass_' + str(i[j]) + '.10'
    Input = str(i[j]) + '.10'
    #tage, r, m, r_star, l_star = np.loadtxt(Input, unpack=True)
    r, r_star, deltar, m, l_star = np.loadtxt(Input, unpack=True)
#    pl.loglog(tage, r/r_sun, label=i[j])
    teff= (l_star/(4*pi*r*r*sigma_b))**0.25
    pl.loglog(m, teff, label=str(i[j]))
    if str(i[j]) == 'nakano_31':
        rteff = (l_star/(4*pi*r_star*r_star*sigma_b))**0.25
        pl.loglog(m, rteff, label=i[j] + 'r')


#pl.xlabel(r'$t_{age}$ $({\rm yr})$', fontsize=18)
pl.xlabel(r'$m$ $({\rm })$', fontsize=18)
pl.ylabel(r'$T_eff$', fontsize=18)
pl.xlim( 4.e32, 2.e36) 
pl.ylim( 1.e3, 1.e5) 
pl.xticks(fontsize=14)
pl.yticks(fontsize=14)
#pl.gcf().subplots_adjust(bottom=0.15)
#pl.rc('text', usetex=False)
#pl.legend(loc="best") 
pl.legend(loc="best", fontsize=14, frameon=False, ncol=2)
pl.savefig('m_vs_teff.pdf')
