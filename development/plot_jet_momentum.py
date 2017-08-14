import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

t16, v16 = np.loadtxt('jet_momentum_16k.txt', unpack=True)
t32, v32 = np.loadtxt('jet_momentum_32k.txt', unpack=True)

plt.rc("text", usetex=True)
plt.rc('font', family='serif')
rcParams.update({'figure.autolayout': True})
rcParams['axes.linewidth']    = 2
rcParams['lines.linewidth']   = 2
rcParams['xtick.major.width'] = 3    # major tick width in points                                                                             
rcParams['xtick.major.size']  = 14    # major tick width in points                                                                            
rcParams['xtick.minor.size']  = 8    # major tick width in points                                                                             
rcParams['xtick.minor.width'] = 2    # major tick width in points                                                                             
rcParams['ytick.major.width'] = 3    # major tick width in points                                                                             
rcParams['ytick.minor.width'] = 2    # major tick width in points                                                                             
rcParams['ytick.major.size']  = 14    # major tick width in points                                                                            
rcParams['ytick.minor.size']  = 8    # major tick width in points                                                                             
rcParams['xtick.labelsize']   = 25
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
plt.minorticks_on()
plt.tick_params('both',length=8, width=1, which='minor')
plt.tick_params('both',length=10, width=1.5, which='major')
plt.tight_layout()

plt.ylabel(r'${\langle v_{\rm jet} \rangle}$ $({\rm km\, s}^{-1})$', fontsize = 25)
plt.xlabel(r'$t-t_*$ $({\rm yr})$', fontsize=25)

plt.loglog(t16, v16, 'b', ls='dashed')
plt.loglog(t32, v32, 'g', label='$u_r$')

plt.xlim(4e4,1e6)
plt.ylim(1e1,1e2)
plt.savefig('jet_v_vs_t.pdf')

