# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

from ctypes import *
from gecatsim.pyfiles.CommonTools import *

def C_set_random_seed(cfg, seed1=None, seed2=None):
    fun = cfg.clib.setall
    fun.argtypes = [c_int, c_int]
    fun.restype = None
    fun(seed1, seed2)
