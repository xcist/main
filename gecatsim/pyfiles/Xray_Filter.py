# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import numpy as np
from scipy import interpolate
from gecatsim.pyfiles.GetMu import GetMu
from gecatsim.pyfiles.CommonTools import *

def Xray_Filter(cfg):
    '''
    Calculate the transmittance of flat and bowtie filters
    return cfg.src.filterTrans, dim [totalNumCells, Ebin]
    Mingye Wu, GE Research
    
    '''
    cfg.src.filterTrans = np.ones([cfg.det.totalNumCells, cfg.spec.nEbin], dtype=np.single)
    cfg = flat_filter(cfg)
    cfg = bowtie_filter(cfg)
    
    return cfg

def flat_filter(cfg):
    '''
    Apply the transmittance of flat filter at source side, dim: [Ebin, totalNumCells]
    '''
    cosineFactors = 1/np.cos(cfg.det.gammas)/np.cos(cfg.det.alphas)
    
    Evec = cfg.sim.Evec
    trans = np.ones([cfg.det.totalNumCells, cfg.spec.nEbin], dtype = np.single)
    
    if hasattr(cfg.protocol, "flatFilter"):
        for ii in range(0, round(len(cfg.protocol.flatFilter)/2)):
            material = cfg.protocol.flatFilter[2*ii]
            depth = cfg.protocol.flatFilter[2*ii+1]
            mu = GetMu(material, Evec)
            trans *= np.exp(-depth*0.1*cosineFactors @ mu.reshape(1, mu.size))
    cfg.src.filterTrans *= trans
    
    return cfg

def bowtie_filter(cfg):
    '''
    Apply the transmittance of bowtie filter, dim: [Ebin, totalNumCells]
    '''
    if not cfg.protocol.bowtie:
        return cfg
    
    # find bowtie file
    bowtieFile = my_path.find("bowtie", cfg.protocol.bowtie, ".txt")

    # read bowtie file
    data = np.loadtxt(bowtieFile, dtype=np.single, comments=['#', '%'])
    
    gammas0 = data[:, 0]
    t0 = data[:, 1:] # thickness in cm
    bowtieMaterials = ['Al', 'graphite', 'Cu', 'Ti']
    
    Evec = cfg.spec.Evec
    gammas1 = cfg.det.gammas
    
    muT = 0
    for i in range(len(bowtieMaterials)):
        mu = GetMu(bowtieMaterials[i], Evec)
        f = interpolate.interp1d(gammas0, t0[:, i], kind='linear', fill_value='extrapolate')
        t1 = f(gammas1)/np.cos(cfg.det.alphas)
        muT += t1 @ mu.reshape(1, mu.size)
    
    trans = np.exp(-muT)
    cfg.src.filterTrans *= trans
    
    return cfg

# if __name__ == "__main__":
#     import matplotlib.pyplot as plt
#     cfg = source_cfg("./cfg/default.cfg")
#
#     cfg = feval(cfg.scanner.detectorCallback, cfg)
#     cfg = feval(cfg.scanner.focalspotCallback, cfg)
#     cfg = feval(cfg.physics.rayAngleCallback, cfg)
#     cfg = feval(cfg.protocol.spectrumCallback, cfg)
#
#     cfg.protocol.bowtie = 'medium'
#     cfg.protocol.flatFilter = ['al', 0.1, 'water', 2]
#
#     cfg = Xray_Filter(cfg)
#     trans = cfg.src.filterTrans.reshape(cfg.scanner.detectorColCount, cfg.scanner.detectorRowCount, cfg.spec.nEbin)
#     check_value(trans)
#     plt.plot(trans[:, 7, 8])
#     plt.show()
    
    
    
    
    
