# Copyright 2025, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import numpy as np
from gecatsim.pyfiles.Spectrum import spectrum_read
from gecatsim.pyfiles.GetMu import GetMu

def CreateHeelEffect(cfg):
    """
    Aim
    Create spectra files that use a simple model for heel effect.
    """
    if not hasattr(cfg, 'heel_effect_limit_angle') or not hasattr(cfg, 'heel_effect_angle_decimation') or not hasattr(
            cfg, 'target_angle'):
        raise ValueError(
            "cfg must have 'heel_effect_limit_angle', 'heel_effect_angle_decimation', and 'target_angle' attributes")

    limitAngle = cfg.heel_effect_limit_angle
    angleSeparation = limitAngle * 2 / (cfg.heel_effect_angle_decimation - 1)

    if limitAngle > cfg.target_angle:
        raise ValueError("cfg.heel_effect_limit_angle should be smaller than cfg.target_angle "
                         "(because the anode(target) can't send x-rays to all angles)")

    x = np.arange(-limitAngle, limitAngle + angleSeparation, angleSeparation)
    nangles = len(x)
    IvecOrig, EvecOrig, TAngleOrig = spectrum_read(cfg.spectrum_filename)

    cfg.spectrum_filename = cfg.spectrum_filename.replace('.dat', '')
    if '\\' in cfg.spectrum_filename:
        folder_char = '\\'
        tmpVar = cfg.spectrum_filename.rfind('\\')
        cfg.spectrum_filename = cfg.spectrum_filename[tmpVar + 1:]
    else:
        folder_char = '/'
        tmpVar = cfg.spectrum_filename.rfind('/')
        cfg.spectrum_filename = cfg.spectrum_filename[tmpVar + 1:]

    intens = np.diag(IvecOrig) @ HeelEffectIntensity(x, EvecOrig, cfg)

    for i in range(nangles):
        if x[i] < 0:
            filename = f"{cfg.spectrum_dir}{folder_char}{cfg.spectrum_filename}_{i:03d}_{x[i]:.0f}.dat"
        else:
            filename = f"{cfg.spectrum_dir}{folder_char}{cfg.spectrum_filename}_{i:03d}_{x[i]:.0f}.dat"

        with open(filename, 'w') as fp:
            fp.write(f"{len(EvecOrig)}\n")
            for j in range(len(EvecOrig)):
                fp.write(f"{EvecOrig[j]},{intens[j, i]}\n")
            fp.write(f"{x[i] + cfg.target_angle + cfg.tube_tilt:.2f}")


def HeelEffectIntensity(theta, Energy, cfg):
    mu = GetMu('w', Energy)

    # Adjust the calculation to prevent overflow
    exp_term = np.clip(-cfg.anode_electron_penetration_in_mm * 0.1 * np.cos(cfg.target_angle / 57.296) *
                       mu * (1 / np.sin((cfg.target_angle / 57.296) -
                                        (cfg.reference_spectrum_angle - cfg.tube_tilt) / 57.296)), -700, 700)
    conversion_factor = (np.ones(mu.shape)) / np.exp(exp_term)

    exp_term_y = np.clip(-cfg.anode_electron_penetration_in_mm * 0.1 * np.cos(cfg.target_angle / 57.296) *
                         mu[:, np.newaxis] * (np.ones(theta.shape) / np.sin(
        (cfg.target_angle / 57.296) + theta / 57.296 - cfg.tube_tilt / 57.296)), -700, 700)
    y = np.diag(np.ones(Energy.shape) * conversion_factor) @ np.exp(exp_term_y)

    return y