from yt.mods import *
import numpy
import argparse

plt_prefix = "BB_hdf5_chk_"

parser = argparse.ArgumentParser(description = "start number to end number")

parser.add_argument('start', metavar='N1', type=int)
parser.add_argument('end', metavar='N2', type=int)

args = parser.parse_args()

time_arr = []
gas_mass_arr = []
star_mass_arr = []
outputArray = [[]]
i_arr = []

tstar = None
for i in range(args.start,args.end) :	
    fn_plt = plt_prefix+"{0:04d}".format(i)
    print "load file " + fn_plt
    pf = load(fn_plt)
    time =pf.h.parameters["time"] 
    data = pf.h.all_data()
    if(tstar == None) : 
	tstar = numpy.min(data["particle_creation_time"])
    
    Msun = 1.9885e33
    rho_mean = 3e-22
    L = 4.8e19
    total_mass = rho_mean*L*L*L/Msun
    star_mass = data["ParticleMassMsun"].sum()
    gas_mass = total_mass-star_mass
    
    time_arr.append(time)
    gas_mass_arr.append(gas_mass)
    star_mass_arr.append(star_mass)
    i_arr.append(i)

print "tstar = " + str(tstar)

numpy.savetxt( "test.out", zip(i_arr, time_arr, star_mass_arr, gas_mass_arr), fmt="%12.7e")
