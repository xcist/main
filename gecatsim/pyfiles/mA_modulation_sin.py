# Copyright 2025, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import numpy as np

def mA_modulation_sin(cfg, viewnr):

    print(f'Calculating mA modulation using {__file__}')

    # Get the rotation angle
    viewangle = cfg.start_angle * np.pi / 180.0 + (viewnr - 1) / cfg.views_per_rotation * 2 * np.pi * cfg.rotation_direction

    # Angle dependent mA modulation scale factor
    mA_scale_factor = 1.0 + cfg.mA_modulation * np.sin(viewangle)

    return mA_scale_factor