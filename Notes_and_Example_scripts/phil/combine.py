import numpy as np
import h5py
import glob
import argparse

plt_prefix = "BB_hdf5_plt_cnt_"
part_prefix = "BB_hdf5_part_"

parser = argparse.ArgumentParser(description = "start number to end number")

parser.add_argument('start', metavar='N1', type=int)
parser.add_argument('end', metavar='N2', type=int)

args = parser.parse_args()

for i in range(args.start,args.end) :	


    fn_plt = plt_prefix+"{0:04d}".format(i)
    fn_part = part_prefix+"{0:04d}".format(i)
    file_plt_exist = glob.glob(fn_plt)
    if not file_plt_exist:
        print 'File: "', fn_plt, '" does not exist, moving to the next file'
        continue
    file_part_exist = glob.glob(fn_plt)
    if not file_part_exist:
        print 'File: "', fn_part, '" does not exist, moving to the next file'
        continue

    print "opening " + fn_plt
    f_plt = h5py.File(fn_plt, "r+")
    print "opening " + fn_part
    f_part = h5py.File(fn_part, "r")

    # localnp
    if(f_part["localnp"][:].sum() == 0) :
         del f_plt, f_part
         continue
    print "processing: {0:04d}".format(i)
     
    f_plt.require_dataset("localnp", f_part["localnp"].shape,
                          f_part["localnp"].dtype)
    
    f_plt["localnp"][:] = f_part["localnp"][:]

    # particle names
        
    f_plt.require_dataset("particle names", f_part["particle names"].shape,
                          f_part["particle names"].dtype)
    
    f_plt["particle names"][:,:] = f_part["particle names"][:,:]
    
    # tracer particles
    
    f_plt.require_dataset("tracer particles", f_part["tracer particles"].shape,
                          f_part["tracer particles"].dtype)
    
    f_plt["tracer particles"][:,:] = f_part["tracer particles"][:,:]
    
    f_plt.close()
    f_part.close()

    del f_plt, f_part
