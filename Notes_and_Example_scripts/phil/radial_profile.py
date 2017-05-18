from yt.mods import * 
import argparse

filename = "BB_hdf5_plt_cnt_0077"
fileout="test.out"
radius = 2.0
pf = load(filename)
sp = pf.h.sphere("max", (radius, "pc"))

p = BinnedProfile1D(sp, 40, "Radiuspc", 1e-3, 3.0, log_space=True)
p.add_fields("Density", weight="CellVolume")
p.add_fields("RadialVelocityKMS", weight="CellMassMsun")
p.add_fields("TangentialOverVelocityMagnitude", weight="CellMassMsun")
p.add_fields("VelocityMagnitude", weight="CellMassMsun")
p.write_out(fileout)
