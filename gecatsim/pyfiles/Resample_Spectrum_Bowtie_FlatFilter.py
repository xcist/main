# Copyright 2025, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import numpy as np

def Resample_Spectrum_Bowtie_FlatFilter(cfg, bowtie, spec, FiltrationTransVec):

    print(f'Resampling "bowtie" and "spec" structures, and filtration transmission...')

    dims = np.array([pad3(bowtie.transVec.shape), pad3(spec.Ivec.shape)])

    if not ((dims[0, 1] == dims[1, 1]) and np.all(bowtie.sinalphas == spec.sinalphas)):  # different fan angles
        print('Resampling "bowtie" and "spec" in the fan (alpha) direction')
        if np.all(dims[:, 1] > 1):
            raise ValueError('We do not yet have code to deal with different fan angle grids. Perhaps now would be a good time to write it.')
        elif dims[0, 1] > 1:  # bowtie depends on fan, need to repmat the spec.Ivec
            spec.Ivec = np.tile(spec.Ivec, (1, dims[0, 1], 1) if len(spec.Ivec.shape) > 2 else (1, dims[0, 1]))
            spec.sinalphas = bowtie.sinalphas
            spec.n_alphas = len(spec.sinalphas)
        elif dims[1, 1] > 1:
            bowtie.transVec = np.tile(bowtie.transVec, (1, dims[1, 1], 1))
            bowtie.sinalphas = spec.sinalphas

    if not ((dims[0, 2] == dims[1, 2]) and np.all(bowtie.singammas == spec.singammas)):  # different cone angles
        print('Resampling "bowtie" and "spec" in the cone (gamma) direction')
        if np.all(dims[:, 2] > 1):
            raise ValueError('We do not yet have code to deal with different cone angle grids. Perhaps now would be a good time to write it.')
        elif dims[0, 2] > 1:  # bowtie depends on cone, need to repmat the spec.Ivec
            spec.Ivec = np.tile(spec.Ivec, (1, 1, dims[0, 2]))
            spec.singammas = bowtie.singammas
            spec.n_gammas = len(spec.singammas)
        elif dims[1, 2] > 1:
            bowtie.transVec = np.tile(bowtie.transVec, (1, 1, dims[1, 2]))
            bowtie.singammas = spec.singammas

    # Compute and store TotalSpectrum: the product of bowtie and filter attenuation and the spectrum
    # (all now a function of energy bin and detector pixel location).
    # The bowtie and filter attenuations are dimensionless.
    # The spectrum is in units of photons/view/mm^2 per energy bin at 1m.
    spec.netIvec = spec.Ivec * bowtie.transVec * FiltrationTransVec

    print('... done resampling.')

    return bowtie, spec, FiltrationTransVec

def pad3(vec):
    vec_pad = np.ones(3)
    vec_pad[:len(vec)] = vec
    return vec_pad