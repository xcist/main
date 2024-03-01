# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

from ctypes import *
from numpy.ctypeslib import ndpointer
from gecatsim.pyfiles.GetMu import GetMu
from gecatsim.pyfiles.CommonTools import *

def C_Phantom_Polygonal_Clear(cfg, num_polygons=None):
    print('Clearing the POLYGONAL phantom in C global variables.')
    func = cfg.clib.clear_polygonalized_phantom
    num_polygons = np.int32(num_polygons)
    func.argtypes = [c_int]
    func.restype = None
    func(num_polygons)