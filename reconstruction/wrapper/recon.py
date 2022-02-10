# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

import numpy as np
import matplotlib.pyplot as plt
import FDK_equiAngle as fdk
from catsim.CommonTools import *


def FDK(cfg):
    # data and parameters
    prep = load_prep(cfg)
    check_value(prep)
    # FDK equianglar
    img = fdk.recon(cfg, prep)
    
    # save image
    
def Helical(cfg):
    pass


    
    
def load_prep(cfg):
    prep = rawread(cfg.resultsName+'.prep', [cfg.protocol.viewCount, cfg.scanner.detectorRowCount, cfg.scanner.detectorColCount], 'float')
    return prep
