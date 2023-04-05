# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

import numpy as np
from catsim.pyfiles.CommonTools import *
#from catsim.pyfiles.Detection_Lag import Detection_Lag

def Detection_DAS(viewIn, viewId, cfg):

    viewOut = viewIn*cfg.scanner.detectionGain
    
    eNoise = np.float32(np.random.randn(viewIn.size)*cfg.sim.eNoise)
    viewOut += eNoise.reshape(viewIn.shape)
    
    return viewOut
