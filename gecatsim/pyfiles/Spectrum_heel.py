# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import os
import glob
import numpy.matlib as nm
from gecatsim.pyfiles.CommonTools import *

def Spectrum_heel(cfg):
    '''
    Read spectral information from a file and return a structure specifying the spectrum.
    The returned spectrum is resampled to new Evec and is scaled to mA and view time.
    The order is Ebin->row->, the dim is python is [col, row, Ebin], i.e. [pixel, Ebin]
    
    Mingye Wu, GE Research
    
    '''
    
    viewTime = cfg.protocol.rotationTime/cfg.protocol.viewsPerRotation/cfg.sim.subViewCount*cfg.protocol.dutyRatio
    specScale = cfg.protocol.spectrumScaling*cfg.protocol.mA*viewTime
    
    if not cfg.spec:
    #if not hasattr(cfg, 'spec'):
        cfg.spec = CFG()
    
    ###------------- monochromatic
    if cfg.physics.monochromatic>0:
        cfg.spec.nEbin = 1
        cfg.spec.Evec = np.array(cfg.physics.monochromatic, dtype=np.single)
        cfg.spec.Ivec = 1.2e6*np.ones([cfg.det.totalNumCells, 1], dtype=np.single)*specScale
        cfg.sim.Evec = cfg.spec.Evec
        return cfg
    
    ###------------- polychromatic
    # Rescale to match CatSim's unit
    # Spectrum is in unit: photons/sec/<area>/<current> at 1-m distance
    # CatSim uses mm and mA, but some spectrum files use cm or A
    if not cfg.protocol.spectrumUnit_mm:
        specScale *= 1e-2
    if not cfg.protocol.spectrumUnit_mA:
        specScale *= 1e-3
    
    # Read spectrum file
    #cfg.protocol.spectrumFilename = my_path.find("spectrum", cfg.protocol.spectrumFilename, "")
    allfiles = glob.glob(cfg.protocol.spectrumFilename+"/*dat")
    Ivec_all, Avec_all = [], []
    # TODO: need to consider the case when not all rows are used
    for thisfile in allfiles:
        Evec0, Ivec0, takeOffAngle = spectrum_read(thisfile)
        Ivec0 *= specScale
    
        Emin = Evec0[0]-0.5*(Evec0[1]-Evec0[0])
        Emax = Evec0[-1]+0.5*(Evec0[-1]-Evec0[-2])
        nEbin = cfg.physics.energyCount
        Ebin = (Emax-Emin)/nEbin
        
        # resample to new Ebin vector
        Evec1 = np.linspace(Emin+Ebin/2, Emax-Ebin/2, nEbin)
        Ivec1 = overlap(Evec0, Ivec0, Evec1)
        inputEnergy = Evec0 @ Ivec0
        outputEnergy = Evec1 @ Ivec1
        Ivec1 *= inputEnergy/outputEnergy
        
        # change dtype
        Evec1 = Evec1.astype(np.float32)
        Ivec1 = Ivec1.astype(np.float32)

        Ivec_all.append(nm.repmat(Ivec1, cfg.scanner.detectorColCount, 1))
        Avec_all.append(takeOffAngle)

    # note that the spectrum should be indexed from large photons (no heel effects, 9 deg) to small photons(heel effects, 5deg)
    # i.e., indexed from cathode to anode
    # sort based on angle from smalle to large
    Ivec_all = [x for _, x in sorted(zip(Avec_all, Ivec_all))]
    Ivec_all = Ivec_all[::-1] # now from large to small to be consistent with gain factor
    # repeat to all pixels
    # based on Detion_Flux, flux dim is [col, row, EBin]
    Ivec_all = np.moveaxis(Ivec_all,1,0)
    Ivec_all = np.vstack(Ivec_all)

    cfg.spec.nEbin = nEbin
    cfg.spec.Evec = Evec1
    cfg.spec.Ivec = Ivec_all
    cfg.sim.Evec = cfg.spec.Evec

    return cfg


def spectrum_read(spectrumFile):
    '''
    Read the spectrum file.
    The output Evec and Ivec are 2-D numpy array with shape [nEbin, 1]
    
    '''
    d0 = []
    for line in open(spectrumFile, 'r'):
        line = line[:-1] # remove the end '\n'
        if line and line[0].isdigit():
            d0.append(line)
            
    nEbin = int(d0[0])
    Evec = []
    Ivec = []
    for ii in range(1, nEbin+1):
        tmp = [float(x.strip()) for x in d0[ii].split(',')]
        Evec.append(tmp[0])
        Ivec.append(tmp[1])

    Evec = np.array(Evec, dtype='single')
    Ivec = np.array(Ivec, dtype='single')
    
    if len(d0)>nEbin+1:
        takeOffAngle = float(d0[nEbin+1])
    else:
        takeOffAngle = 0
        
    return Evec, Ivec, takeOffAngle


# if __name__ == "__main__":
#
#     cfg = source_cfg("./cfg/default.cfg")
#
#     cfg.det.totalNumCells = 5;
#
#     cfg.protocol.spectrumFilename = "tungsten_tar7.0_120_filt.dat"
#     cfg.physics.energyCount = 10
#     cfg.protocol.spectrumScaling = 1
#     cfg.physics.monochromatic = -1
#
#     cfg.protocol.mA = 200
#     cfg.protocol.rotationTime = 1
#     cfg.protocol.viewsPerRotation = 984
#
#     cfg.sim.subViewCount = 1
#     cfg.protocol.dutyRatio = 1
#
#     #Evec, Ivec, takeOffAngle = spectrum_read(cfg.protocol.spectrumFilename)
#     #check_value(Evec)
#     #check_value(Ivec)
#     #check_value(takeOffAngle)
#
#     cfg = Spectrum(cfg)
#     check_value(cfg.spec.nEbin)
#     check_value(cfg.spec.Evec)
#     check_value(cfg.spec.Ivec)
    
