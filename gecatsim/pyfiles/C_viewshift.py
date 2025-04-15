# Copyright 2025, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

from gecatsim.pyfiles.CommonTools import *

def C_viewshift(cfg,rows, cols, planes, pmax, view, ndx, coef, output):
    '''
        Aim
        Wrapper for call to c function by the same or similar name
    '''
    fun = cfg.clib.viewshift
    fun.restype = ctypes.c_int
    fun.argtypes = [
        ctypes.c_int,  # rows
        ctypes.c_int,  # cols
        ctypes.c_int,  # planes
        ctypes.c_int,  # pmax
        ctypes.c_double,  # view
        ctypes.POINTER(ctypes.c_int),  # ndx
        ctypes.POINTER(ctypes.c_double),  # coef
        ctypes.POINTER(ctypes.c_double)  # output
    ]

    # Convert the ndx and coef lists to ctypes arrays
    ndx_array = (ctypes.c_int * len(ndx))(*ndx)
    coef_array = (ctypes.c_double * len(coef))(*coef)
    output_array = (ctypes.c_double * len(output))(*output)

    # Call the C function viewshift
    result = fun(rows, cols, planes, pmax, view, ndx_array, coef_array, output_array)

    if result != 0:
        raise Exception("Error occurred while calling C function viewshift.")

    return list(output_array)
