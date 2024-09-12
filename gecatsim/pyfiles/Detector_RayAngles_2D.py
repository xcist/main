# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import numpy.matlib as nm
from gecatsim.pyfiles.CommonTools import *

def Detector_RayAngles_2D(cfg):
    '''
    Calculates the ray angles (alpha, beta, gamma) for all the cells of 2D detector.
    alpha: cone angle
    beta:  incident angle
    gamma: fan angle
    rayDistance: sdd per cell
    The order of total cells is row->col, in python the shape is [col, row]. The dim is 2D, [pixel, 1]
    
    Mingye Wu, GE Research
    
    '''
    
    det = cfg.det # shortcut
    xyzSrc = cfg.src.samples.mean(axis=0)
    
    tanGammas = np.zeros([det.totalNumCells, 1], dtype=np.single)
    tanAlphas = np.zeros([det.totalNumCells, 1], dtype=np.single)
    betas = np.zeros([det.totalNumCells, 1], dtype=np.single)
    rayDistance = np.zeros([det.totalNumCells, 1], dtype=np.single)
    
    for m in range(0, det.nMod):
        startInd = det.startIndices[m]
        nCells = det.nCells

        xyzDet = nm.repmat(det.modCoords[m, :], nCells, 1) + \
            det.cellCoords @ (np.c_[det.uvecs[m, :], det.vvecs[m, :]].T)
        xyzR = xyzDet - nm.repmat(xyzSrc, nCells, 1);

        cellInd = range(startInd, startInd+nCells)
        rayDistance[cellInd] = vectornorm(xyzR.T)
                
        tanGammas[cellInd] = make_col(xyzR[:, 0]/abs(xyzR[:, 1]))
        tanAlphas[cellInd] = make_col(xyzR[:, 2]/np.sqrt(np.square(xyzR[:, 0])+np.square(xyzR[:, 1])))
        
        wvec = np.cross(det.uvecs[m, :], det.vvecs[m, :])        
        xyzR = xyzR/nm.repmat(vectornorm(xyzR.T), 1, 3)
        betasTmp = np.arccos(np.minimum(xyzR @ wvec, 1))
        
        isLeftMod = xyzR[:, 0]<wvec[0] # betas of the left mod are opposite to the right mod
        betasTmp[isLeftMod] *= -1
        betas[cellInd] = make_col(betasTmp)
    
    cfg.det.rayDistance = rayDistance
    
    #alphas = np.arctan(tanAlphas)
    cfg.det.alphas = np.arctan(tanAlphas)
    #cfg.det.sinAlphas = np.sin(alphas)
    #cfg.det.cosAlphas = np.cos(alphas)
    #cfg.det.tanAlphas = tanAlphas
    
    cfg.det.betas = betas
    #cfg.det.sinBetas = np.sin(betas)
    #cfg.det.cosBetas = np.cos(betas)
    #cfg.det.tanBetas = np.tan(betas)
    
    #gammas = np.arctan(tanGammas)
    cfg.det.gammas = np.arctan(tanGammas)
    #cfg.det.sinGammas = np.sin(gammas)
    #cfg.det.cosGammas = np.cos(gammas)
    #cfg.det.tanGammas = tanGammas
    
    return cfg


# if __name__ == "__main__":
#     import matplotlib.pyplot as plt
#     cfg = source_cfg("./cfg/default.cfg")
#
#     cfg = feval("Detector_ThirdgenCurved", cfg)
#     cfg = feval("Source_Uniform", cfg)
#     cfg = Detector_RayAngles_2D(cfg)
#
#     check_value(cfg.det.rayDistance)
#     check_value(cfg.det.alphas)
#     check_value(cfg.det.cosBetas)
#     check_value(cfg.det.gammas)
#
#     ga = cfg.det.gammas.reshape(cfg.scanner.detectorColCount, cfg.scanner.detectorRowCount)
#     plt.plot(ga[:, 7])
#     plt.show()
