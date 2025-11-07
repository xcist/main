# Copyright 2025, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

# -----------------------------------------------------------------------
# Aim
#   This code can be used to find a rotation matrix R such that
#
#     R*x1=y1
#     R*x2=y2
#
#   assuming that:
#     1) x1, y1, x2, y2 are all unit vectors.
#     2) the angle between x1 and x2 is the same as the angle between y1 and y2
#
#   The code normalizes all vectors so that condition 1 is satisfied
#   If condition 2 is violated, we solve the orthogonal Procrustes problem
#   to make sure that R*x1 and R*x2 are as close as possible to y1 and y2 respectively.
# -----------------------------------------------------------------------

import numpy as np

def find_R(x1, x2, y1, y2):
    x1 = x1 / np.linalg.norm(x1)
    x2 = x2 / np.linalg.norm(x2)
    x3 = np.cross(x1, x2)

    y1 = y1 / np.linalg.norm(y1)
    y2 = y2 / np.linalg.norm(y2)
    y3 = np.cross(y1, y2)

    A = np.dot(np.array([x1, x2, x3]).T, np.array([y1, y2, y3]))
    U, S, Vt = np.linalg.svd(A)
    R = np.dot(Vt.T, U.T)
    return R