# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license
"""
Aim
    Calculate u oversampling coordinates and steps at the given number of oversamples
"""
import os
from gecatsim.pyfiles.CommonTools import *

def USampling(cfg, total_n, intensive_n):
    colsize = cfg.scanner.detectorColSize
    intensive_len = cfg.scanner.colIntensiveOversampleLength

    # Account depth-dependent magnification of the collimator
    # Assume each collimator in front of scintillator is aligned with x-ray source
    if hasattr(cfg.scanner, 'plateThickness') and cfg.scanner.plateThickness > 0:
        platethickness = cfg.scanner.plateThickness
        plateheight = cfg.scanner.plateHeight
        plateairgap = cfg.scanner.plateAirgap
        shadow_width = platethickness / (cfg.scanner.sdd - plateheight - plateairgap) * cfg.scanner.sdd
    else:
        shadow_width = 0
        cfg.scanner.plateLocation = 0

    # Active region in u depends on collimator location and sizes of collimator and kerf
    if (hasattr(cfg.scanner, 'detectionOnKerf') and cfg.scanner.detectionOnKerf == 1) or cfg.scanner.callbackDetector == 'Detector_SVCT':
        uactive_collimator_edge = colsize - shadow_width
        uactive_collimator_center = colsize - shadow_width
    else:
        uactive_collimator_edge = colsize - max(shadow_width, cfg.scanner.colCast)
        uactive_collimator_center = colsize - shadow_width - cfg.scanner.colCast

    if intensive_n == 0 or intensive_len == 0:
        n_u = total_n
        if cfg.scanner.plateLocation == 0:  # Edge collimator
            uactive = uactive_collimator_edge
            du = uactive / n_u
            us = ((np.arange(1, n_u + 1) - (n_u + 1) / 2) * du)
        elif cfg.scanner.plateLocation == 1:  # Center collimator
            uactive = uactive_collimator_center
            du = uactive / n_u
            if not n_u % 2:
                us = ((np.arange(1, n_u + 1) - (n_u + 1) / 2) * du)
                us = us + np.sign(us) * shadow_width / 2
            elif total_n == 1:
                us = np.array([0])
            else:
                raise ValueError('Detector oversample should be EVEN when the plate is at the center')
        else:
            raise ValueError('Plate has to be either at the center or at the edge')

        steplen = np.ones(n_u) * du
    else:
        du_i = intensive_len / intensive_n
        nu_i = intensive_n
        nu_s = total_n - 2 * nu_i
        if nu_s < 1 or colsize - max(shadow_width, cfg.scanner.colCast) < 2 * intensive_len:
            raise ValueError('Intensive oversampling number/area is bigger than total')
        else:
            if cfg.scanner.plateLocation == 0:  # Edge collimator
                uactive = uactive_collimator_edge
                du_s = (uactive - 2 * intensive_len) / nu_s
                us = np.zeros(total_n)
                us[0:nu_i] = ((np.arange(1, nu_i + 1) - 0.5) * du_i - uactive / 2)
                us[total_n - nu_i:total_n] = us[0:nu_i] + uactive - intensive_len
                us[nu_i:total_n - nu_i] = ((np.arange(1, nu_s + 1) - 0.5) * du_s - (uactive - 2 * intensive_len) / 2)
            elif cfg.scanner.plateLocation == 1:  # Center collimator
                if not nu_s % 2:
                    uactive = uactive_collimator_center
                    du_s = (uactive - 2 * intensive_len) / nu_s
                    us = np.zeros(total_n)
                    us[0:nu_i] = ((np.arange(1, nu_i + 1) - 0.5) * du_i - uactive / 2)
                    us[total_n - nu_i:total_n] = us[0:nu_i] + uactive - intensive_len
                    us[nu_i:total_n - nu_i] = ((np.arange(1, nu_s + 1) - 0.5) * du_s - (uactive - 2 * intensive_len) / 2)
                    us = us + np.sign(us) * shadow_width / 2
                elif total_n == 1:
                    us = np.array([0])
                else:
                    raise ValueError('Detector oversample should be EVEN when the plate is at the center')
            else:
                raise ValueError('Plate has to be either at the center or at the edge')
        steplen = np.zeros(total_n)
        steplen[0:nu_i] = np.ones(nu_i) * du_i
        steplen[total_n - nu_i:total_n] = steplen[0:nu_i]
        steplen[nu_i:total_n - nu_i] = np.ones(nu_s) * du_s

    return us, steplen