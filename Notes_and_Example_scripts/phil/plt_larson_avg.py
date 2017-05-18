import argparse
import matplotlib.pyplot as p
import numpy

parser = argparse.ArgumentParser(description = "file index")
parser.add_argument('fileindex', metavar='N', type=int, nargs='+')
args = parser.parse_args()

vtnacc = 0
vtacc = 0
vmagacc = 0
vmagnacc = 0
vrmsacc = 0
rbin = 0
rbin2 = 0
mtacc = 0
vracc = 0

initialized = False
numTrials = 0

for fileindex in args.fileindex : 
	arr = numpy.loadtxt("larson_{0:02d}.out".format(fileindex))
	rbin = arr[:,0]
	vrmsbin = arr[:,1]
	
	arr = numpy.loadtxt("rad_profile_{0:02d}.out".format(fileindex))
	vtnbin = arr[:,3]
	vtbin = arr[:,2]
	vrbin = arr[:,1]
	rbin2 = arr[:,0]
	vmagbin = arr[:,5]
	vmagnbin = arr[:,6]
	mtbin = arr[:,7]

	numTrials = numTrials + 1

	if(initialized) :
		vtnacc = vtnacc + vtnbin
		vtacc = vtacc + vtbin
		vmagnacc = vmagnacc + vmagnbin
		vmagacc = vmagacc + vmagbin
		vrmsacc = vrmsacc + vrmsbin
		mtacc = mtacc + mtbin - mtbin[0]
		vracc = vracc + vrbin
	else :
		vtnacc = vtnbin
		vtacc = vtbin
		vmagnacc = vmagnbin
		vmagacc = vmagbin
		vrmsacc = vrmsbin
		mtacc = mtbin - mtbin[0]
		vracc = vrbin
		initialized = True
 
p.loglog(rbin,vrmsacc/numTrials/1e5, label="Larson rms")
p.loglog(rbin2,vtacc/numTrials/1e5, label="$v_{rms}(r)$ mass-weighted", linestyle='dashed')
p.loglog(rbin2,vtnacc/numTrials/1e5, label="$v_{rms}(r)$ vol-weighted", linestyle='dashed')
p.loglog(rbin2,vmagnacc/numTrials/1e5, label="$v_{mag}(r)$ vol-weighted", linestyle=':')
p.loglog(rbin2,vmagacc/numTrials/1e5, label="$v_{mag}(r)$ mag-weighted", linestyle=':')
p.loglog(rbin,rbin**0.5, label="$v \\propto d^{1/2}$")
p.loglog(rbin2, -vracc/numTrials/1e5, label="$v_r$", linestyle="-.")
#p.loglog(rbin2, mtacc/numTrials, label="$M(r)$", linestyle="-")
#p.loglog(rbin2, 1e2*rbin2, label="$M(r)\\propto r$", linestyle="-.")
p.xlim(5e-2,6e0)
p.ylim(1e-2,1e1)
p.xlabel("r,d [pc]")
p.ylabel("v [km/s]")
#p.legend(loc=3)
p.savefig("larson_avg.pdf")
