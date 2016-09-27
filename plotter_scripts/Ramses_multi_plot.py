import matplotlib
matplotlib.use("Agg")
from matplotlib import rc
rc("text", usetex=True)
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
rcParams['axes.linewidth'] = 2
rcParams['lines.linewidth'] = 2
rcParams['xtick.major.width'] = 3    # major tick width in points
rcParams['xtick.major.size'] = 14    # major tick width in points
rcParams['xtick.minor.size'] = 8    # major tick width in points
rcParams['xtick.minor.width'] = 2    # major tick width in points
rcParams['ytick.major.width'] = 3    # major tick width in points
rcParams['ytick.minor.width'] = 2    # major tick width in points
rcParams['ytick.major.size'] = 14    # major tick width in points
rcParams['ytick.minor.size'] = 8    # major tick width in points
rcParams['xtick.labelsize']=25

import matplotlib.pyplot as pl
import numpy as np
import matplotlib.backends.backend_pdf as mat_pdf
import argparse
import glob
import sys
import os
import math
import h5py
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA

mp = 1.6e-24
pi = np.pi
parsec = 3e18
sec_in_yr = np.pi* 1e7
Msun = 2e33
G = 6.67e-8
c_s = 2.64e4

#sink_04_filename = 'ramses_jeans_04/sink_00132.info'
#ID_04 = np.loadtxt( sinkfile, usecols=[0], unpack=True, skiprows=3, comments="=")
#ID_04 = np.loadtxt( sink_04_filename, usecols=[0], unpack=True, skiprows=3, comments="=")
#ID_04
#sink_08_filename = 'ramses_jeans_08/sink_00138.info'
#ID_08 = np.loadtxt( sink_08_filename, usecols=[0], unpack=True, skiprows=3, comments="=")
#ID_08
#sink_16_filename = 'ramses_jeans_16/sink_00148.info'
#ID_16 = np.loadtxt( sink_16_filename, usecols=[0], unpack=True, skiprows=3, comments="=")
#ID_16
#sink_32_filename = 'ramses_jeans_32/sink_00141.info'
#ID_32 = np.loadtxt( sink_32_filename, usecols=[0], unpack=True, skiprows=3, comments="=")
#sink_32_filename = 'ramses_jeans_32/sink_00139.info'
#ID_32 = np.loadtxt( sink_32_filename, usecols=[0], unpack=True, skiprows=3, comments="=")
#ID_32

Text_04_filename = 'ramses_jeans_04/rad_profile_0132_shellsphere_1.0.out'
Text_08_filename = 'ramses_jeans_08/rad_profile_0138_shellsphere_1.0.out'
Text_16_filename = 'ramses_jeans_16/rad_profile_0144_shellsphere_1.0.out' # 144 is 2.7 # 148 is 4 solar mass
Text_32_filename = 'ramses_jeans_32/rad_profile_0139_shellsphere_1.0.out'
Text_64_filename = 'ramses_jeans_32c/rad_profile_0105_shellsphere_1.out'




rbin64, vrbin64, vrmsbin64, vrmsnbin64, vKbin64, vmagbin64, vmagnbin64, mTbin64, rhobin64, mdotbin64, norm64, angXbin64, angYbin64, angZbin64, vphi_magbin64, sum_of_velocities64, partMass64, part_creation_time64, current_time64, vrms_r_bin64, vrms_l_bin64, vrms_theta_bin64, vrms_phi_bin64 = np.loadtxt(Text_64_filename, unpack=True)
rbin4, vrbin4, vrmsbin4, vrmsnbin4, vKbin4, vmagbin4, vmagnbin4, mTbin4, rhobin4, mdotbin4, norm4, angXbin4, angYbin4, angZbin4, vphi_magbin4, sum_of_velocities4, partMass4, part_creation_time4, current_time4, vrms_r_bin4, vrms_l_bin4, vrms_theta_bin4, vrms_phi_bin4 = np.loadtxt(Text_04_filename, unpack=True)
rbin8, vrbin8, vrmsbin8, vrmsnbin8, vKbin8, vmagbin8, vmagnbin8, mTbin8, rhobin8, mdotbin8, norm8, angXbin8, angYbin8, angZbin8, vphi_magbin8, sum_of_velocities8, partMass8, part_creation_time8, current_time8, vrms_r_bin8, vrms_l_bin8, vrms_theta_bin8, vrms_phi_bin8 = np.loadtxt(Text_08_filename, unpack=True)
rbin16, vrbin16, vrmsbin16, vrmsnbin16, vKbin16, vmagbin16, vmagnbin16, mTbin16, rhobin16, mdotbin16, norm16, angXbin16, angYbin16, angZbin16, vphi_magbin16, sum_of_velocities16, partMass16, part_creation_time16, current_time16, vrms_r_bin16, vrms_l_bin16, vrms_theta_bin16, vrms_phi_bin16 = np.loadtxt(Text_16_filename, unpack=True)
rbin32, vrbin32, vrmsbin32, vrmsnbin32, vKbin32, vmagbin32, vmagnbin32, mTbin32, rhobin32, mdotbin32, norm32, angXbin32, angYbin32, angZbin32, vphi_magbin32, sum_of_velocities32, partMass32, part_creation_time32, current_time32, vrms_r_bin32, vrms_l_bin32, vrms_theta_bin32, vrms_phi_bin32 = np.loadtxt(Text_32_filename, unpack=True)


pl.clf()
pl.rc('text', usetex=True)
pl.rc('font', family='serif')
pl.loglog( rbin4, vrmsbin4/1e5, 'g.-', label='$N_{J}=4$')
pl.loglog( rbin8, vrmsbin8/1e5, 'b', label='$N_{J}=8$')
pl.loglog( rbin16, vrmsbin16/1e5, 'c--', label='$N_{J}=16$')
pl.loglog( rbin32, vrmsbin32/1e5, 'r^-', label='$N_{J}=32$')

#pl.loglog( rbin64, -4.0*pi*rhobin64*vrbin64*(parsec*rbin64)**2*sec_in_yr/Msun, label='$v_T105$')
#pl.loglog( rbin4, -4.0*pi*rhobin4*vrbin4*(parsec*rbin4)**2*sec_in_yr/Msun, 'g.-', label='$v_T171$')
#pl.loglog( rbin8, -4.0*pi*rhobin8*vrbin8*(parsec*rbin8)**2*sec_in_yr/Msun, label='$v_T150$')
#pl.loglog( rbin16, -4.0*pi*rhobin16*vrbin16*(parsec*rbin16)**2*sec_in_yr/Msun, label='$v_T145$')
#pl.loglog( rbin32, -4.0*pi*rhobin32*vrbin32*(parsec*rbin32)**2*sec_in_yr/Msun, label='$v_T131$')
pl.axhline( 2.64e4/1e5, color = 'black',label='$c_s$')
pl.ylim(1e-1, 3e1)
pl.xlim(3e-3, 3e0)
pl.xticks(fontsize=24)
pl.yticks(fontsize=24)
pl.gcf().subplots_adjust(bottom=0.15)
pl.ylabel(r'$v_T$ $({\rm km\, s}^{-1})$', fontsize = 25)
pl.xlabel(r'$r$ $({\rm pc})$', fontsize=25)

pl.rc('text', usetex=False)
pl.legend(loc=1, fontsize=21, frameon=False, ncol=2)
pl.savefig('Ramses_velocity_conversion.pdf')
