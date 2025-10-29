# Copyright 2025, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

# -----------------------------------------------------------------------
# Aim
#   Returns the isotropic 2D Transfer Function of Lorentzian distribution
#   a / (a + f^2)
# -----------------------------------------------------------------------

import numpy as np

def GetLorentzian(ptch, kernelsize, cfg):
    H = np.zeros((kernelsize, kernelsize))
    Ny = kernelsize // 2 + 1  # Nyquist frequency - Positive spectral limit

    a = cfg.lorentzian_a
    b = cfg.lorentzian_b
    c = cfg.lorentzian_c

    u = np.concatenate((np.arange(0, kernelsize // 2 + 1), np.flip(np.arange(1, kernelsize // 2)))) / kernelsize / ptch
    v = np.concatenate((np.arange(0, kernelsize // 2 + 1), np.flip(np.arange(1, kernelsize // 2)))) / kernelsize / ptch

    # Rotation --------------------------------------------------
    r2 = np.tile(u, (kernelsize, 1)) ** 2 + (np.tile(v, (kernelsize, 1)) ** 2).T
    H = (c * a / (a + r2) + (1 - c) * b / (b + r2)) ** 1

    return H