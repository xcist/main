# Copyright 2025, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

# -----------------------------------------------------------------------
# Aim
#   Return an image/volume of radial distances to a central point/axis
# -----------------------------------------------------------------------
import numpy as np

def rimg(nrcols, nrrows=None, nrplanes=1, centercol=None, centerrow=None, centerplane=None):
    if nrrows is None:
        nrrows = nrcols
    if centercol is None:
        centercol = (nrcols + 1) / 2.0
    if centerrow is None:
        centerrow = (nrrows + 1) / 2.0
    if centerplane is None:
        centerplane = (nrplanes + 1) / 2.0

    cols = np.arange(1, nrcols + 1) - centercol
    col2img = np.ones((nrrows, 1)) * (cols ** 2)
    rows = np.arange(1, nrrows + 1).reshape(-1, 1) - centerrow
    row2img = (rows ** 2) * np.ones((1, nrcols))
    radius2img = col2img + row2img
    rimg = np.sqrt(radius2img)
    return rimg