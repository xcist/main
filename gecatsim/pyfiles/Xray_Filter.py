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
    Apply the transmittance of flat filter at source side, optimized version
    '''
    cosineFactors = 1 / (np.cos(cfg.det.gammas) * np.cos(cfg.det.alphas))
    
    Evec = cfg.sim.Evec
    trans = np.ones([cfg.det.totalNumCells, cfg.spec.nEbin], dtype=np.single)
    
    if hasattr(cfg.protocol, "flatFilter"):
        # Pre-allocate array for mu values
        mu_total = np.zeros_like(Evec, dtype=np.single)
        
        for ii in range(0, len(cfg.protocol.flatFilter) // 2):
            material = cfg.protocol.flatFilter[2*ii]
            depth = cfg.protocol.flatFilter[2*ii+1]
            mu = GetMu(material, Evec)
            mu_total += depth * 0.1 * mu
        
        # Single exponential calculation instead of multiple
        trans *= np.exp(-np.outer(cosineFactors, mu_total))
    
    cfg.src.filterTrans *= trans
    return cfg

def bowtie_filter(cfg):
    '''
    Apply the transmittance of bowtie filter, optimized version
    '''
    if not cfg.protocol.bowtie:
        return cfg
    
    bowtieFile = my_path.find("bowtie", cfg.protocol.bowtie, ".txt")
    data = np.loadtxt(bowtieFile, dtype=np.single, comments=['#', '%'])
    
    gammas0 = data[:, 0]
    t0 = data[:, 1:]
    bowtieMaterials = ['Al', 'graphite', 'Cu', 'Ti']
    
    Evec = cfg.spec.Evec
    gammas1 = cfg.det.gammas
    cos_alphas_inv = 1.0 / np.cos(cfg.det.alphas)
    
    # Pre-compute all mu values
    mu_array = np.array([GetMu(material, Evec) for material in bowtieMaterials], dtype=np.single)
    
    # Vectorized interpolation for all materials at once
    f_interps = [interpolate.interp1d(gammas0, t0[:, i], kind='linear', fill_value='extrapolate', assume_sorted=True) 
                 for i in range(len(bowtieMaterials))]
    
    # Compute total attenuation in one operation
    muT = np.zeros((len(gammas1), len(Evec)), dtype=np.single)
    for i, f in enumerate(f_interps):
        t1 = f(gammas1) * cos_alphas_inv
        muT += np.outer(t1, mu_array[i])
    
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
    
    
    
    
    
