import numpy as np
import matplotlib.pyplot as pl

from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA

mp = 1.6e-24
Msun = 2e33
# rhobin is in number density
# avg_density_hd8_hd8.out
rbin_hd4, rhobin_hd4 = np.loadtxt( 'avg_density_hd4_0109.out', usecols=[0,1], unpack=True, comments="=")
rbin_hd8, rhobin_hd8 = np.loadtxt( 'avg_density_hd8_0110.out', usecols=[0,1], unpack=True, comments="=")
rbin_hd32, rhobin_hd32 = np.loadtxt( 'avg_density_hd32_0082.out', usecols=[0,1], unpack=True, comments="=")

print rhobin_hd4
print ""
print rhobin_hd8
print ""
print rhobin_hd32
print ""
#print rhobin_hd4/Msun
#
#rhobin_hd4 = rhobin_hd4/Msun
#rhobin_hd8 = rhobin_hd8/Msun
#rhobin_hd32 = rhobin_hd32/Msun
#
pl.rc('text', usetex=True)
pl.rc('font', family='serif')
host = host_subplot(111, axes_class=AA.Axes)
par1 = host.twinx()
Ndensity_Y_min = 1e1
Ndensity_Y_max = 1e8
host.set_xlim(3e-3, 4e0)
host.set_ylim(Ndensity_Y_min, Ndensity_Y_max)
Mdensity_Y_min = Ndensity_Y_min * 2. * mp
Mdensity_Y_max = Ndensity_Y_max * 2. * mp
par1.set_ylim(Mdensity_Y_min, Mdensity_Y_max)
par1.set_yscale('log')
pl.xticks(fontsize=30)
pl.yticks(fontsize=30)

pl.gcf().subplots_adjust(bottom=0.15)
host.set_ylabel('$n$ $({\\rm cm^{-3}})$', fontsize=25)
host.set_xlabel('$r$ $({\\rm pc})$', fontsize=25)
par1.set_ylabel('$\\rho$ $({\\rm g \\, cm}^{-3})$', fontsize=25)
host.axis["left"].label.set_fontsize(25)
host.axis["bottom"].label.set_fontsize(25)
par1.axis["right"].label.set_fontsize(25)
# rhobin should be in g/cm**3
#host.loglog(rbin, rhobin/mp, label='Number density') 	

host.loglog(rbin_hd4, rhobin_hd4, 'r', label=r'N_J = 4')
host.loglog(rbin_hd8, rhobin_hd8, 'g', label=r'N_J = 8')
host.loglog(rbin_hd32, rhobin_hd32, 'b', label=r'N_J = 32')

pl.legend(loc='best')
pl.rc('text', usetex=False)
pl.savefig('avg_density_convergence.pdf')
