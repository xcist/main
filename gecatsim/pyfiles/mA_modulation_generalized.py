# Copyright 2025, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import numpy as np
def mA_modulation_generalized(cfg, viewnr):
    print(f'Calculating mA modulation using {__name__}')

    if hasattr(cfg, 'max_mA') and hasattr(cfg, 'min_mA'):
        viewangle = cfg.start_angle * np.pi / 180.0 + (viewnr - 1.0) / cfg.views_per_rotation * 2.0 * np.pi * cfg.rotation_direction
        actual_mA = (cfg.max_mA + cfg.min_mA) / 2.0 + cfg.mA_modulation * (cfg.max_mA - cfg.min_mA) / 2.0 * np.cos(2.0 * viewangle)
    else:
        mod_table = np.load(cfg.modulation_table_file_path)['mod_table']
        iView = np.mod((viewnr - 1), cfg.views_per_rotation) + 1
        iRot = np.floor((viewnr - 1) / cfg.views_per_rotation) + 1
        actual_mA = mod_table[iView - 1, int(iRot) - 1]

    mA_scale_factor = actual_mA / cfg.mA
    return mA_scale_factor