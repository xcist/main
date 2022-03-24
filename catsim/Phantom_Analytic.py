# Generated with SMOP  0.41
#from libsmop import *
# Phantom_Analytic.m

# -----------------------------------------------------------------------
#   Program Name: Phantom_Analytic.m                                       
#   Authors:  Samit Basu, Bruno De Man, Jed Pack (GE Global Research)
#   Organization: GE Global Research and GE Healthcare, General Electric Company
#   Version:  6.0.3
#   Date:  Feb 3, 2015
#   Class: GE Confidential. General Electric Proprietary Data (c) 2012 General Electric Company
    
    # Aim
#   This function is basically a wrapper for ParsePhantom_pp_or_ppm
    
    # Inputs
#   cfg                  The configuration structure
    
    # Outputs
#   Materials            Cell array of phantom materials (strings)
#   NumberOfMaterials    Number of materials in phantom
    
    # History: 
#   2012-10-09 Paul FitzGerald (GE Global Research)
#              Renamed - was GetPhantom.m.
#              Cleaned up and added "Verbose" output. 
# -----------------------------------------------------------------------

import numpy as np
from ctypes import *
from numpy.ctypeslib import ndpointer
import os, json
from catsim.GetMu import GetMu
from catsim.CommonTools import *
from catsim.Phantom_Analytic_Get import Phantom_Analytic_Get
    
def Phantom_Analytic(cfg):
    
    print("Starting to read ANALYTIC phantom...")

    phobj, phobject = Phantom_Analytic_Get(cfg)
    #NumberOfMaterials = len(Materials)

    cfg.phantom.numberOfMaterials = len(phobject['materialList'])
    set_materials(cfg, phobject['materialList'])
    #TODO: need to read c source code to finish this part
    #set_volume(cfg, phobj)
    set_detector(cfg)
    set_source(cfg, phobj)

    #print('Phantom contains {%d} materials'.format(NumberOfMaterials))
    print('... done reading phantom.')
    #return Materials,NumberOfMaterials

    #should return cfg
    return cfg

# TODO: need materialsList passed in here
def set_materials(cfg, materialList):
    Evec = cfg.spec.Evec
    nMat = len(materialList)
    Mus = np.zeros([Evec.size, nMat], dtype=np.float64)
    for i in range(nMat):
        Mus[:, i] = GetMu(materialList[i], Evec)/10

    # source file: base_prop/src/analytic_projector.h
    fun = cfg.clib.set_material_info
    fun.argtypes = [c_int, c_int, ndpointer(c_double)]
    fun.restype = None
    fun(nMat, Evec.size, Mus)
    #breakpoint()

#TODO: not done yet
# what excatly does set_volume mean?
# in C: void set_phantom_info(int numObjs, int *objType, int *clipStInd, int *nPlanes, int *matInd, double *objCent, double *shp, dou     ble *Qmat, double *clipNormVec, double *clipDist, double *den, int totalNumPlanes)
# in matlabe: calllib(CatSimLib, 'set_phantom_info', numObjects, T, cumCP, NumCP, materialIndex, X, K, Q, Eta, S, D, cumCP_end);
def set_volume(cfg, phobj):
    fun = cfg.clib.set_phantom_info
    fun.argtypes = [POINTER(c_int), ndpointer(c_float), ndpointer(c_int), \
        c_float, c_float, c_float, c_float, c_float, ndpointer(c_ubyte), c_int, c_int]
    fun.restype = None
    
    Status = [0]
    Status = (c_int*1)(*Status)
    fun(Status, volumeData, volumeDims, \
        offsets[0], offsets[1], offsets[2], voxelsize[0], voxelsize[2], xyMask, materialIndex, numberOfMaterials)

def set_source(cfg):
    SourceWeights = cfg.src.weights
    NumberOfSubSources = cfg.src.nSamples
    
    fun = cfg.clib.set_src_info_vox
    fun.argtypes = [ndpointer(c_float), c_int]
    fun.restype = None    
    fun(SourceWeights, NumberOfSubSources)

def set_detector(cfg):
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
    fun = cfg.clib.set_module_info
    fun.argtypes = [POINTER(c_float), POINTER(c_float), POINTER(c_int), ndpointer(c_float), \
                POINTER(c_int), ndpointer(c_float), ndpointer(c_float), c_int, c_int, c_int, c_int]
    fun.restype = None
    fun(Height, Width, Pix, Coords, \
                Sub, Sampling, Weight, nModuleTypes, maxPix, maxSubDets, moduleOverlapType)
    
