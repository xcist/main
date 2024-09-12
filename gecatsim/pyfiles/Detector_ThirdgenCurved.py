# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import numpy.matlib as nm
import math
from gecatsim.pyfiles.CommonTools import *

def Detector_ThirdgenCurved(cfg):
    '''
    Returns a modular 3rd generation detector (multi-row) focused at the source.
    The order of total cells is row->col, in python the shape is [col, row].The dim is 2D, [pixel, xyz]
    
    Mingye Wu, GE Research
    
    '''
    
    # shortcuts
    sid = cfg.scanner.sid
    sdd = cfg.scanner.sdd
    nRow = cfg.scanner.detectorRowsPerMod
    nCol = cfg.scanner.detectorColsPerMod
    nMod = math.ceil(cfg.scanner.detectorColCount/nCol)
    rowSize = cfg.scanner.detectorRowSize
    colSize = cfg.scanner.detectorColSize
    
    # cell coords
    cols = (np.arange(0, nCol)-(nCol-1)/2)*colSize
    rows = (np.arange(0, nRow)-(nRow-1)/2)*rowSize
    cols = nm.repmat(cols, nRow, 1).T.reshape(nCol*nRow, 1)
    rows = nm.repmat(rows, 1, nCol).T
    cellCoords = np.c_[cols, rows]
    
    # sample U coords
    nu = cfg.physics.colSampleCount
    du = colSize*cfg.scanner.detectorColFillFraction/nu
    us = (np.arange(0, nu)-(nu-1)/2)*du
        
    # sample V coords
    nv = cfg.physics.rowSampleCount
    dv = rowSize*cfg.scanner.detectorRowFillFraction/nv
    vs = (np.arange(0, nv)-(nv-1)/2)*dv
    
    # sample coords
    nSamples = nu*nv
    us = nm.repmat(us, nv, 1).T.reshape(nSamples, 1)
    vs = nm.repmat(vs, 1, nu).T
    sampleCoords = np.c_[us, vs]
    weights = np.ones((nu, nv))/nSamples
    
    # module offset
    modWidth = nCol*colSize
    dAlpha = 2*math.atan(modWidth/2/sdd)
    uOffset = math.atan(cfg.scanner.detectorColOffset*colSize/sdd)
    vOffset = cfg.scanner.detectorRowOffset*rowSize
    alphas = (np.arange(0, nMod)-(nMod-1.0)/2)*dAlpha + uOffset
    alphas = make_col(alphas)
    
    # module coords, uvecs, vvecs
    sinA = np.sin(alphas)
    cosA = np.cos(alphas)
    modCoords = np.c_[sdd*sinA, sid-sdd*cosA, nm.repmat(vOffset, nMod, 1)]
    uvecs = np.c_[cosA, sinA, np.zeros((nMod, 1))]
    vvecs = np.c_[(np.zeros((nMod, 1)), np.zeros((nMod, 1)), np.ones((nMod, 1)))]
    startIndices = np.arange(0, nMod)*nRow*nCol
    
    # detector definition
    #if not hasattr(cfg, 'det'):
    if not cfg.det:
        cfg.det = CFG()
    
    cfg.det.nCells = nRow*nCol
    cfg.det.cellCoords = np.single(cellCoords)
    
    cfg.det.nSamples = nSamples
    cfg.det.sampleCoords = np.single(sampleCoords)
    cfg.det.weights = np.single(weights)
    cfg.det.activeArea = colSize*cfg.scanner.detectorColFillFraction*rowSize*cfg.scanner.detectorRowFillFraction
        
    cfg.det.nMod = nMod
    cfg.det.modCoords = np.single(modCoords)
    cfg.det.uvecs = np.single(uvecs)
    cfg.det.vvecs = np.single(vvecs)
    
    cfg.det.totalNumCells = cfg.scanner.detectorColCount*cfg.scanner.detectorRowCount
    cfg.det.startIndices = np.int32(startIndices)
    cfg.det.nModDefs = 1
    cfg.det.modTypes = np.zeros((nMod, 1), dtype=np.int32)
    
    cfg.det.width = (nCol+1)*colSize
    cfg.det.height = (nRow+1)*rowSize
        
    return cfg


# if __name__ == "__main__":
#
#     cfg = source_cfg("./cfg/default.cfg")
#
#     cfg.scanner.sid = 540
#     cfg.scanner.sdd = 950
#     cfg.scanner.detectorRowsPerMod = 16
#     cfg.scanner.detectorColsPerMod = 16
#     cfg.scanner.detectorColOffset = -1.25
#     cfg.scanner.detectorRowOffset = 0
#     cfg.scanner.detectorColSize = 1.0
#     cfg.scanner.detectorRowSize = 1.0
#     cfg.scanner.detectorColFillFraction = 0.8
#     cfg.scanner.detectorRowFillFraction = 0.8
#     cfg.scanner.detectorColCount = 900
#     cfg.scanner.detectorRowCount = cfg.scanner.detectorRowsPerMod
#
#     cfg.physics.colSampleCount = 4
#     cfg.physics.rowSampleCount = 3
#
#     cfg = Detector_ThirdgenCurved(cfg)
#     check_value(cfg.det.cellCoords)
#     check_value(cfg.det.sampleCoords)
#     check_value(cfg.det.modCoords)
#     check_value(cfg.det.uvecs)
    
