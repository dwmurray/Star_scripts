import argparse
import matplotlib 
matplotlib.use("Agg")
import numpy
from yt.mods import *
import matplotlib.pyplot as p

parser = argparse.ArgumentParser(description = "file index")
parser.add_argument('start', metavar='N1', type=int)
parser.add_argument('end', metavar='N1', type=int)
args = parser.parse_args()

for i in range(args.start,args.end+1) : 
	pf=load("BB_hdf5_plt_cnt_{0:04d}".format(i))
	dd=pf.h.all_data()

	time = pf.current_time - dd["particle_creation_time"]
	mass = dd["ParticleMassMsun"]
	time = time/(6.7e6*3.15e7)

	p.scatter(time, mass)
	p.scatter(time.max(), mass.sum(), s=50)

t = numpy.arange(0.,0.3,0.01)
p.plot(t, 8000.*t*t, label="$8000 t^2$")
p.legend(loc="best")

p.savefig("mass_part.pdf")
