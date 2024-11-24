# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import os
import ctypes
from gecatsim.pyfiles.CommonTools import *

def load_catsim_lib():
    lib_path = my_path.paths["lib"]
    my_path.add_dir_to_path(lib_path)

    # load C/C++ lib
    ll = ctypes.cdll.LoadLibrary
    if os.name == "nt":
        libFile = "libcatsim64.dll"
    else:
        libFile = "libcatsim.so"
    clib = ll(os.path.join(lib_path, libFile))

    return clib

def set_rndseed(seed):

    if seed < 0:
        print('Warning: Seed is set to negative, will be converted to positive')

    seed = abs(round(seed)) + 1  # needs to be positive integer

    clib = load_catsim_lib()
    clib.setall(ctypes.c_int32(seed * 10), ctypes.c_int32(seed))