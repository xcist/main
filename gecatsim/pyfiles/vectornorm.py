# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license
"""
Aim
    Compute the norms of an array of vectors
"""
import numpy as np

def vectornorm(xyz):
    if xyz.shape[0] != 3:
        print('ERROR: argument of vectornorm has to be of size 3 x n')
        return None

    norms = np.sqrt(np.sum(xyz * xyz, axis=0))
    return norms