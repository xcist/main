
from ctypes import *
from numpy.ctypeslib import ndpointer
from gecatsim.pyfiles.CommonTools import *
import gecatsim.pyfiles.CommonTools as c

# Assuming you have loaded the shared library CatSimLib using ctypes.
# For example:
# CatSimLib = ctypes.CDLL('CatSimLib.dll')  # Replace 'CatSimLib.dll' with the actual library file.

def C_Phantom_Polygonal_SetPolygon(cfg, vertices = None, numTriangles = None, density = None, ID = None):
    print(f'Setting a polygon in C global variables; density = {density}; ID = {ID}.')
    func = cfg.clib.pass_polygon_to_c
    cVertices = vertices.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    # Assuming vertices is a list of floats, and num_triangles and ID are integers.
    func.argtypes = [
        ctypes.POINTER(ctypes.c_double),  # vertices as a double pointer (C array of doubles)
        ctypes.c_int,                     # num_triangles as an integer
        ctypes.c_double,                  # density as a double
        ctypes.c_int                      # ID as an integer
    ]

    # Call the C function and get the return value
    func.restype = None
    func(cVertices, numTriangles, density, ID)

    return cVertices
