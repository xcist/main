# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

import re
import sys
from ctypes import *
from numpy.ctypeslib import ndpointer
from gecatsim.pyfiles.GetMu import GetMu
from gecatsim.pyfiles.CommonTools import *

# python conversion from C_Materials_Polygonal_Set.m
def set_materials(cfg, materialList):
    Evec = np.array(cfg.spec.Evec)
    nMat = len(materialList)
    Mus = np.zeros([Evec.size, nMat], dtype=np.float64)
    for i in range(nMat):
        Mus[:, i] = GetMu(materialList[i], Evec)/10

    # polygonal_projector.c: void set_material_info(int materialCount, int eBinCount, double *muTable)
    fun = cfg.clib.set_material_info
    fun.argtypes = [c_int, c_int, ndpointer(c_double)]
    fun.restype = None
    fun(nMat, Evec.size, Mus)

def C_Phantom_Polygonal_Clear(cfg, num_polygons = None):
    print('Clearing the POLYGONAL phantom in C global variables.')
    func = cfg.clib.clear_polygonalized_phantom
    num_polygons= np.int32(num_polygons)
    func.argtypes = [c_int]
    func.restype = None
    func(num_polygons)