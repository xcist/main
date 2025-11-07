# Copyright 2025, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import numpy as np
from gecatsim.pyfiles.GetMu import GetMu
from scipy.interpolate import splev, splrep

def HVLDoseCalculations(cfg, spec):
    """
    Calculate central HVL and central air kerma of spectrum
    """
    Evec = spec.Evec[:, 0]
    fromNist = [4.742, 1.334, 5.389e-1, 1.537e-1, 6.833e-2, 4.098e-2, 3.041e-2, 2.407e-2, 2.325e-2, 2.496e-2, 2.672e-2]
    fromNistX = [10, 15, 20, 30, 40, 50, 60, 80, 100, 150, 200]
    mutrbyro = splev(Evec, splrep(fromNistX, fromNist))  # in cm^2/g
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
        col += 1
        row += 1
        netIvec = (netIvec + spec.netIvec[:, row + cfg.row_count * (col - 1)]) / 2

    if hasattr(cfg, 'calculate_central_dose') and cfg.calculate_central_dose:
        if not hasattr(cfg, 'dose_distance'):
            cfg.dose_distance = 1000

        flux = netIvec * np.exp(-air_mu * t_air * 0.1) * (1000 / r) ** 2  # in Photons/mm^2 at r
        tmp = flux * mutrbyro * Evec * 100 * 1000 * 1.60218e-16 / 87 / 1e-4 * 1e3  # 100mm^2/cm^2, 1000g/kg, 1.60218e-16J/keV, 87*10^-4Gr/R, 1000mR/R

        dose = np.sum(tmp)
        print(f'\nbCT Air kerma at {r} mm = {dose} mR ({dose * 87 * 1e-4} mGy) , after {t_air} mm of air')

        r = cfg.dose_distance
        t_air = 0
        flux = netIvec * np.exp(-air_mu * t_air * 0.1) * (1000 / r) ** 2  # in Photons/mm^2 at r
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
            dose = 0
            flux = netIvec * np.exp(-air_mu * (r - thick[j]) * 0.1 - al_mu * thick[j] * 0.1) * (
                        1000 / r) ** 2  # in Photons/mm^2 at r
            tmp = flux * mutrbyro * Evec * 100 * 1000 * 1.60218e-16 / 87 / 1e-4 * 1e3
            dose = np.sum(tmp)
            D[j] = dose

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


