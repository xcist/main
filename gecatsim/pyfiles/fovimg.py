# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import numpy as np

def fovimg(nrrows, nrcols, nrplanes=1, radius=None, centerrow=None, centercol=None):
    # Set default input arguments if not provided
    if radius is None:
        radius = min(nrrows, nrcols) / 2.0
    if centerrow is None:
        centerrow = (1 + nrrows) / 2.0
    if centercol is None:
        centercol = (1 + nrcols) / 2.0

    # Calculate squared row and column indices; then find where radius is smaller than input radius
    row2img = (np.arange(1, nrrows + 1)[:, None] - centerrow) ** 2 * np.ones((1, nrcols))
    col2img = np.ones((nrrows, 1)) * (np.arange(1, nrcols + 1) - centercol) ** 2

    img = (row2img + col2img) <= radius ** 2

    # Replicate nrplanes times if 3D volume (cylinder) is requested
    if nrplanes > 1:
        img = np.repeat(img[:, :, np.newaxis], nrplanes, axis=2)

    return img
