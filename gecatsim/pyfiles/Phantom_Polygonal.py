# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

import numpy as np
import os
from gecatsim.pyfiles.Phantom_Analytic import parse_analytical_ppm
from gecatsim.pyfiles.Phantom_Polygonal_ReadPolygon import extract_polygonal_objects
from gecatsim.pyfiles.CommonTools import *
from ctypes import *
from numpy.ctypeslib import ndpointer
from gecatsim.pyfiles.GetMu import GetMu
from gecatsim.pyfiles.CommonTools import *

def Phantom_Polygonal(cfg):
    print('Starting to pass POLYGONAL phantom to C...')

    cfg.phantom.filename = my_path.find('phantom', cfg.phantom.filename, '')
    # if os.path.exists(cfg.phantom.filename):
    #    exec(open(cfg.phantom.filename).read())

    object = extract_polygonal_objects(cfg.phantom.filename)

    numPoly = len(object['type'])
    feval('C_Phantom_Polygonal_Clear',cfg,numPoly)

    for index in range(numPoly):
        # rotate around y-axis
        if hasattr(cfg.phantom, "phantomRotationY") and cfg.phantom.phantomRotationY != 0:
            ang = cfg['phantom_rotation_y']
            sina = np.sin(-ang / 180 * np.pi)
            cosa = np.cos(-ang / 180 * np.pi)
            transmat = np.array([[cosa, 0, sina], [0, 1, 0], [-sina, 0, cosa]])
            object['vertices'][index] = np.dot(object['vertices'][index], transmat)

        # Scale and position
        object['vertices'][index] = object['vertices'][index] * cfg.phantom.scale + np.tile(cfg.phantom.centerOffset, (object['vertices'][index].shape[0], 1))
        # Pass phantom info to C
        num_triangles = object['num_triangles'][index]
        _triangles = []
        for vertex_list in object['vertices'][index]:
            for i in range(0,len(vertex_list),3):
                triangle = np.array(vertex_list[i:i+3])
                _triangles.extend(triangle.flatten())
        triangles = np.array(_triangles).reshape((9, num_triangles))

        # triangles = object['vertices'][index][inds.ravel(), :].T.reshape((9, num_triangles))
        feval('C_Phantom_Polygonal_SetPolygon', cfg, triangles, triangles.shape[1], object['density'],
              object['materialId'])

    print('... done with phantom.')
    return cfg

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

