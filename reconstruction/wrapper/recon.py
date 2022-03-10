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
    img = feval(cfg.scanner.recontype, cfg, prep)

    # save image
    save_raw(cfg, img)

    
def Helical(cfg):
    pass

    
def load_prep(cfg):
    prep = rawread(cfg.resultsName+'.prep', [cfg.protocol.viewCount, cfg.scanner.detectorRowCount, cfg.scanner.detectorColCount], 'float')
    return prep


def save_raw(cfg, matrix):

    if cfg.recon.unit =='hu':
        matrix = matrix*(1000/(cfg.recon.mu))
        matrix = matrix + cfg.recon.hu_offset
    elif cfg.recon.unit == 'mm':
        matrix = matrix
    elif cfg.recon.unit == 'cm':
        matrix = matrix*10

    cfg.resultsBasename = cfg.phantom.filename
    fname = cfg.resultsBasename+'.recon'
    # rawwrite(fname, matrix)
    matrix.tofile(fname + '.raw')