# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

from ctypes import *
from numpy.ctypeslib import ndpointer
from gecatsim.pyfiles.GetMu import GetMu
from gecatsim.pyfiles.CommonTools import *

#
# Pass det and src to C projectors
#
def C_Projector_SetData(cfg, viewId):
    if not cfg.sim.isPhantomScan:
        return
    
    if viewId == cfg.sim.startViewId or cfg.physics.recalcDet:
        projectorIDs = get_projector_id(cfg)
        for projectorID in projectorIDs:
            func = 'set_detector_' + projectorID
            eval(func)(cfg)
            
    if viewId == cfg.sim.startViewId or cfg.physics.recalcSrc:
        projectorIDs = get_projector_id(cfg)
        for projectorID in projectorIDs:
            func = 'set_source_' + projectorID
            eval(func)(cfg)

def get_projector_id(cfg):
    if not isinstance(cfg.phantom.callback, list):
        phantom_callbacks = [cfg.phantom.callback]
    else:
        phantom_callbacks = cfg.phantom.callback
    projectorIDs = []
    for theCallback in phantom_callbacks:
        if 'analytic' in theCallback.lower():
            projectorIDs.append('analytic')
        elif 'voxelized' in theCallback.lower():
            projectorIDs.append('voxelized')
        elif any(substring in theCallback.lower() for substring in ['ncat','polygon']):
            projectorIDs.append('ncat')
    return projectorIDs

#
# Analytic Projector
#
def set_detector_analytic(cfg):
    det = cfg.det
    
    Height = [det.height]
    Width = [det.width]
    Pix = [det.nCells]
    Coords = det.cellCoords.astype(np.double)
    Sub = [det.nSamples]
    Sampling = det.sampleCoords.astype(np.double)
    Weight = det.weights.astype(np.double)
    nModuleTypes = det.nModDefs
    maxPix = np.max(det.nCells)
    maxSubDets = np.max(det.nSamples)
    moduleOverlapType = 2
    
    Height = (c_double*1)(*Height)
    Width = (c_double*1)(*Width)
    Pix = (c_int*1)(*Pix)
    Sub = (c_int*1)(*Sub)
        
    # analytic_projector.c: void set_module_info(double *Height, double *Width, int *Pix, double *Coords, int *Sub, double *Sampling, double *Weight, int nModuleTypes, int maxPix, int maxSubDets, int moduleOverlapType)
    fun = cfg.clib.set_module_info
    fun.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_int), ndpointer(c_double), \
        POINTER(c_int), ndpointer(c_double), ndpointer(c_double), c_int, c_int, c_int, c_int]
    fun.restype = None
    fun(Height, Width, Pix, Coords, \
        Sub, Sampling, Weight, nModuleTypes, maxPix, maxSubDets, moduleOverlapType)

def set_source_analytic(cfg):
    # analytic_projector.c: void set_src_info(double *sourceWeights, int nSubSources)
    fun = cfg.clib.set_src_info
    fun.argtypes = [ndpointer(c_double), c_int]
    fun.restype = None
    
    SourceWeights = cfg.src.weights.astype(np.double)
    NumberOfSubSources = cfg.src.nSamples
    
    fun(SourceWeights, NumberOfSubSources)

#
# Voxelized Projector
#
def set_detector_voxelized(cfg):
    det = cfg.det
    
    Height = [det.height]
    Width = [det.width]
    Pix = [det.nCells]
    Coords = det.cellCoords
    Sub = [det.nSamples]
    Sampling = det.sampleCoords
    Weight = det.weights
    nModuleTypes = det.nModDefs
    maxPix = np.max(det.nCells)
    maxSubDets = np.max(det.nSamples)
    moduleOverlapType = 2
    
    Height = (c_float*1)(*Height)
    Width = (c_float*1)(*Width)
    Pix = (c_int*1)(*Pix)
    Sub = (c_int*1)(*Sub)
    
    # the C func wants data order: row -> col -> pixel_ind or sample_ind
    fun = cfg.clib.set_module_info_vox
    fun.argtypes = [POINTER(c_float), POINTER(c_float), POINTER(c_int), ndpointer(c_float), \
        POINTER(c_int), ndpointer(c_float), ndpointer(c_float), c_int, c_int, c_int, c_int]
    fun.restype = None    
    fun(Height, Width, Pix, Coords, \
        Sub, Sampling, Weight, nModuleTypes, maxPix, maxSubDets, moduleOverlapType)
 
def set_source_voxelized(cfg):
    SourceWeights = cfg.src.weights
    NumberOfSubSources = cfg.src.nSamples
    
    fun = cfg.clib.set_src_info_vox
    fun.argtypes = [ndpointer(c_float), c_int]
    fun.restype = None    
    fun(SourceWeights, NumberOfSubSources)
 
#
# NCAT and Polygonal Projector
#
def set_detector_ncat(cfg):
    det = cfg.det
    
    Height = [det.height]
    Width = [det.width]
    Pix = [det.nCells]
    Coords = det.cellCoords.astype(np.double)
    Sub = [det.nSamples]
    Sampling = det.sampleCoords.astype(np.double)
    Weight = det.weights.astype(np.double)
    nModuleTypes = det.nModDefs
    maxPix = np.max(det.nCells)
    maxSubDets = np.max(det.nSamples)
    UNUSED = 0
    
    Height = (c_double*1)(*Height)
    Width = (c_double*1)(*Width)
    Pix = (c_int*1)(*Pix)
    Sub = (c_int*1)(*Sub)
        
    # in nCAT_main.c: void set_module_info_NCAT(double *Height, double *Width, int *Pix, double *Coords, int *Sub, double *Sampling, double *Weight, int nModuleTypes, int maxPix, int maxSubDets, int UNUSED)
    # the C func wants data order: row -> col -> pixel_ind or sample_ind
    fun = cfg.clib.set_module_info_NCAT
    fun.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_int), ndpointer(c_double), \
        POINTER(c_int), ndpointer(c_double), ndpointer(c_double), c_int, c_int, c_int, c_int]
    fun.restype = None    
    fun(Height, Width, Pix, Coords, \
        Sub, Sampling, Weight, nModuleTypes, maxPix, maxSubDets, UNUSED)
 
def set_source_ncat(cfg):
    # in nCAT_main.c: void set_src_info_NCAT(double *sourceWeights, int nSubSources)
    fun = cfg.clib.set_src_info_NCAT
    fun.argtypes = [ndpointer(c_double), c_int]
    fun.restype = None
    
    SourceWeights = cfg.src.weights.astype(np.double)
    NumberOfSubSources = cfg.src.nSamples
    
    fun(SourceWeights, NumberOfSubSources)
 
