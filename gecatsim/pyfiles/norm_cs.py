# Copyright 2025, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import numpy as np
# This file is used as a replacement for Freemat's norm function
#  It is needed because Freemat has a bug on 64-bit Linux builds prior to version 3.6

G_fm_norm_is_ok = None

def norm_cs(in_array):
    global G_fm_norm_is_ok

    if G_fm_norm_is_ok is None:
        try:
            ver = np.__version__
        except:
            ver = '2.0'

        if float(ver[:3]) < 3.6:
            G_fm_norm_is_ok = False
        else:
            G_fm_norm_is_ok = True

    if not G_fm_norm_is_ok:
        if in_array.ndim == 1 or (in_array.ndim == 2 and 1 in in_array.shape):
            return np.sqrt(np.sum(in_array * in_array))
        else:
            raise ValueError('norm_cs expects a vector')
    else:
        return np.linalg.norm(in_array)