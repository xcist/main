# Copyright 2025, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

# -----------------------------------------------------------------------
# Aim
#   Calculate central HVL and central air kerma of spectrum
# -----------------------------------------------------------------------

import numpy as np
from gecatsim.pyfiles.GetMu import GetMu
def HVLDoseCalculations(cfg, spec):
    Evec = spec.Evec[:, 0]
    fromNist = [4.742, 1.334, 0.5389, 0.1537, 0.06833, 0.04098, 0.03041, 0.02407, 0.02325, 0.02496, 0.02672]
    fromNistX = [10, 15, 20, 30, 40, 50, 60, 80, 100, 150, 200]
    mutrbyro = np.interp(Evec, fromNistX, fromNist)  # in cm^2/g
    r = cfg.dose_distance
    t_air = r
    al_mu = GetMu('al', Evec)
    air_mu = GetMu('air', Evec)

    if cfg.col_count % 2:
        col = (cfg.col_count + 1) // 2
        row = (cfg.row_count + 1) // 2
        netIvec = spec.netIvec[:, row + cfg.row_count * (col - 1)]
    else:
        col = cfg.col_count // 2
        row = cfg.row_count // 2
        netIvec = spec.netIvec[:, row + cfg.row_count * (col - 1)]
        if col < cfg.col_count and row < cfg.row_count:
            col += 1
            row += 1
            if row + cfg.row_count * (col - 1) < spec.netIvec.shape[1]:
                netIvec = (netIvec + spec.netIvec[:, row + cfg.row_count * (col - 1)]) / 2

    netIvec = netIvec.reshape(-1)  # Ensure netIvec is a 1D array

    # Interpolate air_mu, mutrbyro, and Evec to match the shape of netIvec
    air_mu = np.interp(Evec, Evec[:len(air_mu)], air_mu)
    mutrbyro = np.interp(Evec, Evec[:len(mutrbyro)], mutrbyro)
    Evec = np.interp(Evec, Evec[:len(Evec)], Evec)

    # Ensure all arrays have the same shape
    if not (netIvec.shape == air_mu.shape == mutrbyro.shape == Evec.shape):
        raise ValueError(f"Shape mismatch: netIvec shape {netIvec.shape}, air_mu shape {air_mu.shape}, mutrbyro shape {mutrbyro.shape}, Evec shape {Evec.shape}")

    if hasattr(cfg, 'calculate_central_dose') and cfg.calculate_central_dose:
        if not hasattr(cfg, 'dose_distance'):
            cfg.dose_distance = 1000

        flux = netIvec * np.exp(-air_mu * t_air * 0.1) * (1000 / r)**2  # in Photons/mm^2 at r
        tmp = flux * mutrbyro * Evec * 100 * 1000 * 1.60218e-16 / 87 / 1e-4 * 1e3  # 100mm^2/cm^2, 1000g/kg, 1.60218e-16J/keV, 87e-4Gr/R, 1000mR/R

        dose = np.sum(tmp)
        print(f'\nbCT Air kerma at {r} mm = {dose} mR ({dose * 87e-4} mGy) , after {t_air} mm of air')

        r = cfg.dose_distance
        t_air = 0
        flux = netIvec * np.exp(-air_mu * t_air * 0.1) * (1000 / r)**2  # in Photons/mm^2 at r
        tmp = flux * mutrbyro * Evec * 100 * 1000 * 1.60218e-16 / 87 / 1e-4 * 1e3
        dose2 = np.sum(tmp)
        print(f'bCT Air kerma at {r} mm = {dose2} mR ({dose2 / 114.5} mGy) , after {t_air} mm of air')

        mu_efec = -np.log(dose / dose2) / r / 0.1
        print(f'bCT Effective air mu for spectra = {mu_efec}')

    if hasattr(cfg, 'calculate_hvl') and cfg.calculate_hvl:
        print('\nlinear Interpolation -----------------------------------')
        nrslices = 400
        thick = 0.1 * np.arange(nrslices + 1)
        D = np.zeros(nrslices + 1)

        for j in range(nrslices + 1):
            flux = netIvec * np.exp(-air_mu * (r - thick[j]) * 0.1 - al_mu * thick[j] * 0.1) * (1000 / r)**2  # in Photons/mm^2 at r
            tmp = flux * mutrbyro * Evec * 100 * 1000 * 1.60218e-16 / 87 / 1e-4 * 1e3
            D[j] = np.sum(tmp)

        Dnorm = D / D[0]
        HVL = np.zeros(4)
        HVL[0] = np.interp(-np.log(0.5), -np.log(Dnorm), thick)
        print(f'HVL = {HVL[0]} mm of Al')

        HVL[1] = np.interp(-np.log(0.25), -np.log(Dnorm), thick)
        print(f'1/4 HVL = {HVL[1]} mm of Al')

        HVL[2] = np.interp(-np.log(0.125), -np.log(Dnorm), thick)
        print(f'1/8 HVL = {HVL[2]} mm of Al')

        HVL[3] = np.interp(-np.log(0.1), -np.log(Dnorm), thick)
        print(f'1/10 HVL = {HVL[3]} mm of Al')