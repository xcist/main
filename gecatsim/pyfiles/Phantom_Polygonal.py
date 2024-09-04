# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import numpy as np
import os
from ctypes import *
from numpy.ctypeslib import ndpointer
from gecatsim.pyfiles.GetMu import GetMu
from gecatsim.pyfiles.CommonTools import *
from gecatsim.pyfiles.Phantom_Polygonal_ReadPolygon import extract_polygonal_objects

def Phantom_Polygonal(cfg):
    ###----------- pass material Mu to C
    set_material(cfg)
    ###----------- pass volume info to C
    set_volume(cfg)
    return cfg


def set_volume(cfg):
    # print('Starting to pass POLYGONAL phantom to C...')

    cfg.phantom.filename = my_path.find('phantom', cfg.phantom.filename, '')
    # if os.path.exists(cfg.phantom.filename):
    #    exec(open(cfg.phantom.filename).read())

    object = extract_polygonal_objects(cfg.phantom.filename)

    numPoly = len(object['type'])
    C_Phantom_Polygonal_Clear(cfg, numPoly)

    for index in range(numPoly):
        # rotate around y-axis
        if hasattr(cfg.phantom, "phantomRotationY") and cfg.phantom.phantomRotationY != 0:
            ang = cfg['phantom_rotation_y']
            sina = np.sin(-ang / 180 * np.pi)
            cosa = np.cos(-ang / 180 * np.pi)
            transmat = np.array([[cosa, 0, sina], [0, 1, 0], [-sina, 0, cosa]])
            object['vertices'][index] = np.dot(object['vertices'][index], transmat)

        # Scale and position
        object['vertices'][index] = object['vertices'][index]*cfg.phantom.scale + np.tile(cfg.phantom.centerOffset,(object['vertices'][index].shape[0], 1)).astype(np.float32)

        # # Move objects to the center, for testing puprpose ONLY.
        # for kk in [2]:
            # object['vertices'][index][:,kk] -= np.mean(object['vertices'][index][:,kk])
            
        # Pass volumes to C
        C_Phantom_Polygonal_SetPolygon(cfg, object['vertices'][index], object['num_triangles'][index], object['density'][index], object['materialId'][index])

    # print('... done with phantom.')
    return cfg

def C_Phantom_Polygonal_Clear(cfg, num_polygons = None):
    # print('Clearing the POLYGONAL phantom in C global variables.')
    func = cfg.clib.clear_polygonalized_phantom
    num_polygons = np.int32(num_polygons)
    func.argtypes = [c_int]
    func.restype = None
    func(num_polygons)

def C_Phantom_Polygonal_SetPolygon(cfg, vertices = None, numTriangles = None, density = None, ID = None):
    # print(f'Setting a polygon in C global variables; density = {density}; ID = {ID}.')
    # In nCAT_main.c: void pass_polygon_to_c(float *vertices, int num_triangles, float density, int ID)
    func = cfg.clib.pass_polygon_to_c
    func.argtypes = [ndpointer(c_float), c_int, c_float, c_int]
    func.restype = None
    func(vertices, numTriangles, density, ID)


def set_material(cfg):
    materialList = [
        'ncat_water',
        'ncat_muscle',
        'ncat_lung',
        'ncat_dry_spine',
        'ncat_dry_rib',
        'ncat_adipose',
        'ncat_blood',
        'ncat_heart',
        'ncat_kidney',
        'ncat_liver',
        'ncat_lymph',
        'ncat_pancreas',
        'ncat_intestine',
        'ncat_skull',
        'ncat_cartilage',
        'ncat_brain',
        'ncat_spleen',
        'ncat_iodine_blood',
        'ncat_iron',
        'ncat_pmma',
        'ncat_aluminum',
        'ncat_titanium',
        'ncat_air',
        'ncat_graphite',
        'ncat_lead',
        'ncat_breast_mammary',
        'ncat_skin',
        'ncat_iodine',
        'ncat_eye_lens',
        'ncat_ovary',
        'ncat_red_marrow',
        'ncat_yellow_marrow',
        'ncat_testis',
        'ncat_thyroid',
        'ncat_bladder']
    
    Evec = cfg.spec.Evec
    nMat = len(materialList)
    Mus = np.zeros([Evec.size, nMat], dtype=np.double)
    for i in range(nMat):
        Mus[:, i] = GetMu(materialList[i], Evec)/10    # cm^-1 --> mm^-1

    # in nCAT_main.c: void set_material_info_polygon(int materialCount, int eBinCount, double *muTable)
    # the C func wants data order: materialIndex -> Ebin, so the dim is [Ebin, materialIndex]
    fun = cfg.clib.set_material_info_polygon
    fun.argtypes = [c_int, c_int, ndpointer(c_double)]
    fun.restype = None
    fun(nMat, Evec.size, Mus)
    
    return nMat

