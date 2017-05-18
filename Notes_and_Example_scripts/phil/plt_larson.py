import argparse
import matplotlib.pyplot as p
import numpy

parser = argparse.ArgumentParser(description = "file index")
parser.add_argument('fileindex', metavar='N1', type=int)
args = parser.parse_args()

arr = numpy.loadtxt("larson_{0:02d}.out".format(args.fileindex))
rbin = arr[:,0]
vrmsbin = arr[:,1]

arr = numpy.loadtxt("rad_profile_{0:02d}.out".format(args.fileindex))
vtnbin = arr[:,3]
vtbin = arr[:,2]
vrbin = arr[:,1]
rbin2 = arr[:,0]
mtbin = arr[:,5]

p.loglog(rbin,vrmsbin/1e5, label="Larson rms")
p.loglog(rbin2,vtbin/1e5, label="$v_{rms}(r)$ mass-weighted")
p.loglog(rbin2,vtnbin/1e5, label="$v_{rms}(r)$ vol-weighted", linestyle='dashed')
p.loglog(rbin,rbin**0.5, label="$v \\propto d^{1/2}$")
#p.loglog(rbin2, mtbin, label="mass")
p.xlabel("r,d [pc]")
p.ylabel("v [km/s]")
p.legend(loc="best")
p.savefig("larson.pdf")
