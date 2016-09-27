import h5py
import numpy
import os

def getValues( particleFile, value) :
        names = numpy.char.strip(particleFile["particle names"].value)
        starParticles = particleFile["tracer particles"].value
        ind = (numpy.char.strip(names) == value)[:,0]
        return starParticles[:, ind].flatten()

def File_opener(fileIndex):
    #particleFile = h5py.File("{0}{1:04d}".format(part_prefix,fileIndex))
    particleFile = h5py.File("{0}{1}".format(part_prefix,fileIndex))
    current_time = particleFile['real scalars'][0][1]
    print current_time
    print particleFile.keys()
# get all the star particles first
    if("tracer particles" in particleFile.keys()) :

        partIndices = getValues(particleFile, "tag").astype("int")
        xps = getValues(particleFile, "posx")
        yps = getValues(particleFile, "posy")
        zps = getValues(particleFile, "posz")
        vx = getValues(particleFile, "velx")
        vy = getValues(particleFile, "vely")
        vz = getValues(particleFile, "velz")
        
        mass_part_solar = getValues(particleFile, "mass")/Msun
        partCreateTime = getValues(particleFile, "creation_time")
        print partIndices, mass_part_solar
        #partIndices = all_data["particle_index"]
        #xps = all_data["particle_position_x"]
        #yps = all_data["particle_position_y"]
        #zps = all_data["particle_position_z"]
        #mps = all_data["ParticleMassMsun"]


# Constants #
mp = 1.6e-24
pi = numpy.pi
parsec = 3e18
sec_in_yr = numpy.pi* 1e7
Msun = 2e33
G = 6.67e-8
c_s = 2.64e4

part_prefix = 'BB_hdf5_part_'
fileIndex = '0106'
print "{0}{1}".format(part_prefix,fileIndex)

#os.chdir('./q1_08')

file_exist = glob.glob('{0}_{1}_{2:04d}_{3}*.out'.format(file_prefix, quad, framestep, compare_files))
if not file_exist:
	print 'File: "{0}_{1}_{2:04d}_{3}*.out" does not exist!'.format(file_prefix, quad, framestep, compare_files)
	continue

particle_list = glob.glob('{0}_{1}_{2:04d}_{3}*.out'.format(file_prefix, quad, framestep, compare_files))
print 'Looking at: {0}_{1}_{2:04d}_{3}*.out'.format(file_prefix, quad, framestep, compare_files)
print 'The total number of star particles in this timestep is: ', len(particle_list)
for j in range(len(particle_list)):
	particle_number = particle_list[j]
	particle_number = particle_number[-11:-4]
	print particle_number[:-2], withParticleIDValue
	if int(particle_number[:-2]) == withParticleIDValue or (withAllParticles):
		print "Plotting particle", j, 'of', len(particle_list)
		particle_file_exist =  glob.glob('{0}_{1}_{2:04d}_{3}{4}.out'.format(file_prefix, quad, framestep, compare_files, particle_number))
		if not particle_file_exist:
			print 'File: "{0}_{1}_{2:04d}_{3}{4}.out" does not exist!'.format(file_prefix, quad, framestep, compare_files, particle_number)
			continue
