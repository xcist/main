# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

from numpy.ctypeslib import ndpointer
from gecatsim.pyfiles.CommonTools import *

#----------------- load C lib
clib = load_C_lib()
fun = clib.rndpoi
fun.argtypes = [ndpointer(ctypes.c_float), ctypes.c_int]
fun.restype = None

def randpf(arr):
    # numpy array only
    arr1 = arr[arr>0]
    fun(arr1, arr1.size)
    arr[arr>0] = arr1
    
    return arr

# if __name__ == '__main__':
#     data0 = np.float32(np.random.random([5, 4])*100)
#     data0[1:4, 1:3] = 0
#
#     check_value(data0)
#     data1 = randpf(data0)
#     check_value(data0)
#     check_value(data1)
