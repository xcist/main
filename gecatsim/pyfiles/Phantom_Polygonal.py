import numpy as np
import os
from gecatsim.pyfiles.Phantom_Analytic import parse_analytical_ppm
from gecatsim.pyfiles.CommonTools import feval
from gecatsim.pyfiles.CommonTools import *


def phantom_polygonal(cfg):
    print('Starting to pass POLYGONAL phantom to C...')

    cfg.phantom.filename = my_path.find('phantom', cfg.phantom.filename, '')
    # if os.path.exists(cfg.phantom.filename):
    #    exec(open(cfg.phantom.filename).read())

    object = parse_analytical_ppm(cfg.phantom.filename)

    numPoly = len(object['type'])
    feval('C_Phantom_Polygonal_Clear', numPoly)

    for index in range(numPoly):
        # rotate around y-axis
        if hasattr(cfg.phantom, "phantomRotationY") and cfg.phantom.phantomRotationY != 0:
            ang = cfg['phantom_rotation_y']
            sina = np.sin(-ang / 180 * np.pi)
            cosa = np.cos(-ang / 180 * np.pi)
            transmat = np.array([[cosa, 0, sina], [0, 1, 0], [-sina, 0, cosa]])
            object['vertices'][index] = np.dot(object['vertices'][index], transmat)

        # Scale and position
        object['vertices'][index] = object['vertices'][index] * cfg['phantom_scale'] + np.tile(
            cfg['phantom_position_offset'], (object['vertices'][index].shape[0], 1))

        # Pass phantom info to C
        num_triangles = object['tri_inds'][index].shape[0]
        inds = object['tri_inds'][index].T
        triangles = object['vertices'][index][inds.ravel(), :].T.reshape((9, num_triangles))
        feval('C_Phantom_Polygonal_SetPolygon', triangles, triangles.shape[1], object['density'][index],
              object['material'][index] - 1)

    print('... done with phantom.')


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
        Mus[:, i] = GetMu(materialList[i], Evec) / 10

    # polygonal_projector.c: void set_material_info(int materialCount, int eBinCount, double *muTable)
    fun = cfg.clib.set_material_info
    fun.argtypes = [c_int, c_int, ndpointer(c_double)]
    fun.restype = None
    fun(nMat, Evec.size, Mus)


def C_Phantom_Polygonal_Clear(cfg, num_polygons=None):
    print('Clearing the POLYGONAL phantom in C global variables.')
    func = cfg.clib.clear_polygonalized_phantom
    num_polygons = np.int32(num_polygons)
    func.argtypes = [c_int]
    func.restype = None
    func(num_polygons)
