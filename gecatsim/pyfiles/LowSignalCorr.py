# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

from gecatsim.pyfiles.CommonTools import *
from scipy.ndimage import convolve1d
from ctypes import *
from numpy.ctypeslib import ndpointer

# low signal correction
def LowSignalCorr(cfg, prep):
    print("Applying LSC...", end="")
    inshape = prep.shape
    lsc_prep = prep.reshape(cfg.protocol.viewCount, cfg.scanner.detectorRowCount, cfg.scanner.detectorColCount)
    # find negative values and replace them with conv data
    if np.min(lsc_prep) <= 0:
        neg_mask = lsc_prep <= 0
        kernel = [0.5, 0, 0.5]
        tmpprep = convolve1d(lsc_prep, kernel, axis=2) # perform conv in the col direction
        lsc_prep[neg_mask] = tmpprep[neg_mask]

    # now perform real LSC
    func = cfg.clib.negative_log
    func.argtypes = [c_int, c_int, ndpointer(c_float), ndpointer(c_float), c_int]
    func.restype = None
    lsc_prep = np.float32(lsc_prep)
    for i in range(lsc_prep.shape[0]):
        thisview = np.copy(lsc_prep[i], order='C')
        thisviewout = np.zeros(thisview.shape, dtype=np.float32)
        func(lsc_prep.shape[1], lsc_prep.shape[2], thisview, thisviewout, 1)
        lsc_prep[i] = thisviewout

    print("done.\n")
    return lsc_prep.reshape(inshape)
