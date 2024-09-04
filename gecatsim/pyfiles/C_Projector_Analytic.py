# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

from ctypes import *
from numpy.ctypeslib import ndpointer
from gecatsim.pyfiles.CommonTools import *

def C_Projector_Analytic(cfg, viewId, subViewId):
    ###------- C function and interface
    # analytic_projector.c: void Projector(double *Paras, double subviewWeight, double *thisView, double *sourcePoints, int nSubSources, 
    #    double *srcHullPoints, int nSrcHullPoints, int *firstDetIndex, int nModulesIn, int *modTypeInds, 
    #    double *Up, double *Right, double *Center, int UNUSED_tvLength)
    
    ###------- Arguments
    det = cfg.detNew
    src = cfg.srcNew
    
    Paras = np.zeros(9, dtype=np.double)
    subviewWeight = 1.
    thisView = np.zeros([det.totalNumCells, cfg.spec.nEbin], dtype=np.double) # buffer for C
    sourcePoints = src.samples.astype(np.double)
    nSubSources = src.nSamples
    srcHullPoints = src.corners.astype(np.double)
    nSrcHullPoints = src.nCorners
    firstDetIndex = det.startIndices
    nModulesIn = det.nMod
    modTypeInds = det.modTypes
    Up = det.vvecs.astype(np.double)
    Right = det.uvecs.astype(np.double)
    Center = det.modCoords.astype(np.double)
    UNUSED_tvLength = cfg.spec.nEbin*det.totalNumCells
    UNUSED = 0
    
    numThreads = cfg.phantom.projectorNumThreads
    if os.name == "nt":
         numThreads = 1
         
    if numThreads>1:
        fun = cfg.clib.Projector_threaded
        fun.argtypes = [ndpointer(c_double), c_double, ndpointer(c_double), ndpointer(c_double), c_int, \
            ndpointer(c_double), c_int, ndpointer(c_int), c_int, ndpointer(c_int), \
            ndpointer(c_double), ndpointer(c_double), ndpointer(c_double), c_int, c_int, c_double]
        fun.restype = None
        ###------- Run C function
        fun(Paras, subviewWeight, thisView, sourcePoints, nSubSources, \
            srcHullPoints, nSrcHullPoints, firstDetIndex, nModulesIn, modTypeInds, \
            Up, Right, Center, UNUSED_tvLength, numThreads, UNUSED)
    
    else:
        fun = cfg.clib.Projector
        fun.argtypes = [ndpointer(c_double), c_double, ndpointer(c_double), ndpointer(c_double), c_int, \
            ndpointer(c_double), c_int, ndpointer(c_int), c_int, ndpointer(c_int), \
            ndpointer(c_double), ndpointer(c_double), ndpointer(c_double), c_int]
        fun.restype = None
        ###------- Run C function
        fun(Paras, subviewWeight, thisView, sourcePoints, nSubSources, \
            srcHullPoints, nSrcHullPoints, firstDetIndex, nModulesIn, modTypeInds, \
            Up, Right, Center, UNUSED_tvLength)
    
    ###------- Apply transmittance
    cfg.thisSubView *= thisView

    return cfg
