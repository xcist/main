import numpy as np
from ctypes import *
from numpy.ctypeslib import ndpointer

# make_vol_ncat(float *volume, int Nx, double xoff, double dx, int Ny, double yoff, double dy, int Nz, double zoff, double dz, int oversampling, int UNUSED_num_volumes, int material_volumes)
def C_Volume_Poly_Get(cfg, volume, Nx, xoff, dx, Ny, yoff, dy, Nz, zoff, dz, oversampling, num_volumes, material_volumes): 
    print('Getting voxelized phantom/image material volume(s) from C')
    func = cfg.clib.make_vol_ncat_polygon
    func.argtypes = [ndpointer(c_float),
                    c_int,
                    c_double,
                    c_double,
                    c_int,
                    c_double,
                    c_double,
                    c_int,
                    c_double,
                    c_double,
                    c_int];
    func.restype = None

    func(volume, Nx, xoff, dx, Ny, yoff, dy, Nz, zoff, dz, oversampling, num_volumes, material_volumes)

    return volume
