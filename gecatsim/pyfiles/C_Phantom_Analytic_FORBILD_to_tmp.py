# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

from ctypes import *
import gecatsim as xc
from gecatsim.pyfiles.CommonTools import *

def C_Phantom_Analytic_FORBILD_to_tmp(Scale, ppPhantomFilename, tmpPhantomFilename):

    clib = load_C_lib()
    func = clib.TranslatePhantom_FORBILD_to_tmp
    func.argtypes = [c_double, c_char_p, c_char_p]
    func.restype = None
    func(Scale, ppPhantomFilename.encode('utf-8'), tmpPhantomFilename.encode('utf-8'))
