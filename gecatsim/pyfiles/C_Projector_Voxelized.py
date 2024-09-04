# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

from ctypes import *
from numpy.ctypeslib import ndpointer
from gecatsim.pyfiles.CommonTools import *

def C_Projector_Voxelized(cfg, viewId, subViewId):
    ###------- C function and interface
    fun = cfg.clib.voxelized_projector
    fun.argtypes = [POINTER(c_int), c_float, ndpointer(c_float), ndpointer(c_float), c_int, \
        POINTER(c_float), c_int, POINTER(c_int), \
        c_int, ndpointer(c_int), ndpointer(c_float), ndpointer(c_float), ndpointer(c_float), \
        c_int, c_int, c_int, c_int, c_float, c_int]
    fun.restype = None
    
    ###------- loop of source samples and materials
    det = cfg.detNew
    src = cfg.srcNew
    
    newSrcWeights = np.array([1], dtype=np.single)
    set_source(cfg, newSrcWeights, 1)  # reset source weight and number to 1
    
    trans = 0
    matPVS = np.zeros([det.totalNumCells, cfg.spec.nEbin], dtype=np.single) # buffer for C
    
    for srcId in range(src.nSamples):
        theSourcePoint = src.samples[srcId,:]
        pValueSpectrum = 0
        
        for matId in range(cfg.phantom.numberOfMaterials):
            Status = [0]
            Status = (c_int*1)(*Status)
            unused1 = 1
            matPVS[:] = 0
            sourcePoints = theSourcePoint
            nSubSources = 1
            unused2 = [2]
            unused2 = (c_float*1)(*unused2)
            unused3 = 3
            unused4 = [4]
            unused4 = (c_int*1)(*unused4)
            nModulesIn = det.nMod
            modTypeInds = det.modTypes
            Up = det.vvecs
            Right = det.uvecs
            Center = det.modCoords
            unused5 = 5
            unused6 = 6
            MaterialIndex = matId+1
            MaterialIndexInMemory = MaterialIndex
            unused7 = 7
            freeTheMemory = 0
            
            if viewId==cfg.sim.stopViewId and subViewId==cfg.sim.subViewCount-1 \
                and srcId==src.nSamples-1 and matId==cfg.phantom.numberOfMaterials-1:
                freeTheMemory = 1
            
            # the C func wants data order: (row -> col) or xyz or Ebin -> pixel_ind or sample_ind
            fun(Status, unused1, matPVS, sourcePoints, nSubSources, \
                unused2, unused3, unused4, \
                nModulesIn, modTypeInds, Up, Right, Center, \
                unused5, unused6, MaterialIndex, MaterialIndexInMemory, unused7, freeTheMemory)
            pValueSpectrum += matPVS
        
        trans += src.weights[0,srcId]*np.exp(-pValueSpectrum)
    
    cfg.thisSubView *= trans

    return cfg

def set_source(cfg, SourceWeights, NumberOfSubSources):
    fun = cfg.clib.set_src_info_vox
    fun.argtypes = [ndpointer(c_float), c_int]
    fun.restype = None    
    fun(SourceWeights, NumberOfSubSources)
