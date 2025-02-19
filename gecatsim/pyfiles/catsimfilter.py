# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license
"""
Aim
    Generate filtered spectra

Input
    cfg filename or cfg structure
        cfg.spectrum_SpectrumFilename
        cfg.filter
    spectra data, .dat files

Output
    filtered spectrum files, form: '\\filtered\\*_filt.dat'
"""
import os
import numpy as np
from gecatsim.pyfiles.Spectrum import spectrum_read
from gecatsim.pyfiles.GetMu import GetMu

def catsimfilter(cfg):

    if len(cfg.flat_filters) == 0:
        raise ValueError('No filter specified.')

    if not hasattr(cfg, 'spectrum_filtered_string'):
        cfg.spectrum_filtered_string = ''

    original_suffix = '_orig'
    filtered_suffix = '_filt' + cfg.spectrum_filtered_string

    # Get spectrum path/name
    spectrum_filename = cfg.spectrum_filename
    if not os.path.exists(spectrum_filename):
        full_spectrum_filename = os.path.join('gecatsim/spectrum', spectrum_filename)
        if not os.path.exists(full_spectrum_filename):
            raise FileNotFoundError(f"*** ERROR: {spectrum_filename} cannot be found!")
        spectrum_filename = full_spectrum_filename

    # Read spectrum
    ivec = []
    evec = []
    if os.path.isdir(spectrum_filename):
        spectrum_filenames = [os.path.join(spectrum_filename, f) for f in os.listdir(spectrum_filename) if
                              f.endswith('.dat')]
        num_spectrum_files = len(spectrum_filenames)
        for file_index in range(num_spectrum_files):
            spectrum_filename_with_path = spectrum_filenames[file_index]
            e, i, _ = spectrum_read(spectrum_filename_with_path)
            ivec.append(i)
            evec.append(e)
    elif os.path.isfile(spectrum_filename):
        num_spectrum_files = 1
        e, i, _ = spectrum_read(spectrum_filename)
        ivec.append(i)
        evec.append(e)

    # Calculate and apply filter transmission
    cfg.sim.Evec = evec[0]  # Assuming all spectra have the same energy vector
    cosineFactors = 1 / np.cos(cfg.det.gammas) / np.cos(cfg.det.alphas)
    trans = np.ones([cfg.det.totalNumCells, cfg.spec.nEbin], dtype=np.single)

    if hasattr(cfg.protocol, "flatFilter"):
        for ii in range(0, round(len(cfg.protocol.flatFilter) / 2)):
            material = cfg.protocol.flatFilter[2 * ii]
            depth = cfg.protocol.flatFilter[2 * ii + 1]
            mu = GetMu(material, cfg.sim.Evec)
            trans *= np.exp(-depth * 0.1 * cosineFactors[:, np.newaxis] @ mu.reshape(1, -1))

    cfg.src.filterTrans = trans

    for file_index in range(num_spectrum_files):
        ivec[file_index] = ivec[file_index] * np.exp(-cfg.src.filterTrans)

    # Write spectrum
    if os.path.isdir(spectrum_filename):
        path_filtered = os.path.join(spectrum_filename, 'filtered')
    elif os.path.isfile(spectrum_filename):
        pathstr = os.path.dirname(spectrum_filename)
        path_filtered = os.path.join(pathstr, 'filtered')

    try:
        os.makedirs(path_filtered)
    except FileExistsError:
        pass

    for file_index in range(num_spectrum_files):
        number_ebins = len(evec[file_index])

        if os.path.isdir(spectrum_filename):
            spectrum_filename_original = os.path.basename(spectrum_filenames[file_index])
        elif os.path.isfile(spectrum_filename):
            spectrum_filename_original = os.path.basename(spectrum_filename)

        spectrum_filename_filtered = spectrum_filename_original.replace(original_suffix, filtered_suffix)

        path_name = os.path.join(path_filtered, spectrum_filename_filtered)
        with open(path_name, 'w') as fid:
            fid.write(f"{number_ebins}\n")
            for energy_index in range(number_ebins):
                fid.write(f"{evec[file_index][energy_index]},{ivec[file_index][energy_index]}\n")