# Copyright 2025, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

def mA_modulation_Rev2Rot(cfg, viewnr):

    print(f'Calculating mA modulation using {__file__}')

    if viewnr <= cfg.views_per_rotation:
        mA_scale_factor = 1
    else:
        mA_scale_factor = 100 / cfg.mA

    return mA_scale_factor