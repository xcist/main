# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

from numpy.ctypeslib import ndpointer
from gecatsim.pyfiles.CommonTools import *

# Function to convert FORBILD phantom
def C_Phantom_Analytic_FORBILD_to_tmp(cfg,scale, pp_phantom_filename, tmp_phantom_filename):
    # Print headline (assuming a similar function exists in Python)
    print("\nConverting a FORBILD phantom in C.\n")

    fun = cfg.clib.TranslatePhantom_FORBILD_to_tmp
    fun.argtypes = [
        ctypes.c_double(scale),
        ctypes.c_char_p(pp_phantom_filename.encode('utf-8')),
        ctypes.c_char_p(tmp_phantom_filename.encode('utf-8'))
    ]
    fun.restype = None
    fun(scale, pp_phantom_filename, tmp_phantom_filename)

    return