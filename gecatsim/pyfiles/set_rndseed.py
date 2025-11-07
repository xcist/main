# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import os
import ctypes
from ctypes import *
from gecatsim.pyfiles.CommonTools import *

def set_rndseed(cfg,seed):
    if seed < 0:
        print('Warning: Seed is set to negative, will be converted to positive')

    seed = abs(round(seed)) + 1  # needs to be positive integer

    fun = cfg.clib.setall
    fun.argtypes = [c_int, c_int]
    fun.restype = None
    fun(ctypes.c_int32(seed * 10), ctypes.c_int32(seed))