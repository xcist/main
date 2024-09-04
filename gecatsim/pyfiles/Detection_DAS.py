# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

#import numpy as np
from gecatsim.pyfiles.CommonTools import *
import sys
#from catsim.pyfiles.Detection_Lag import Detection_Lag

def Detection_DAS(viewIn, viewId, cfg):

    # to be consisten with spectrum_heel, gain factors should be indexed from cathode end to anode end
    if hasattr(cfg.scanner.detectionGain, "__len__"):
        if len(cfg.scanner.detectionGain) != cfg.scanner.detectorRowCount:
            print("Error! detectionGain should have the same length as row count\n")
            sys.exit(1)
        else:
            # based on plot, viewIn is flatten of [col, row]
            tmp = np.tile(cfg.scanner.detectionGain, cfg.scanner.detectorColCount)
            viewOut = np.float32(viewIn*tmp)
    else:
        viewOut = viewIn*cfg.scanner.detectionGain
    
    eNoise = np.float32(np.random.randn(viewIn.size)*cfg.sim.eNoise)
    viewOut += eNoise.reshape(viewIn.shape)
    
    return viewOut
