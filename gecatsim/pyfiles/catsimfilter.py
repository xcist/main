# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license
import os
import numpy as np
from gecatsim.pyfiles.Spectrum import spectrum_read
from gecatsim.pyfiles.Xray_Filter import flat_filter

def catsimfilter(cfg):

    if len(cfg.flat_filters) == 0:
        raise ValueError('No filter specified.')

    if not hasattr(cfg, 'spectrum_filtered_string'):
        cfg.spectrum_filtered_string = ''

    original_suffix = '_orig'
    filtered_suffix = '_filt' + cfg.spectrum_filtered_string

    # Get spectrum path/name
    spectrum_filename = os.path.exists(cfg.spectrum_filename)
    if not spectrum_filename:
        full_spectrum_filename = os.path.join('gecatsim/spectrum', cfg.spectrum_filename)
        spectrum_filename = os.path.exists(full_spectrum_filename)
        if not spectrum_filename:
            raise FileNotFoundError(f"*** ERROR: {cfg.spectrum_filename} cannot be found!")

    # Read spectrum
    ivec = []
    evec = []
    if os.path.isdir(spectrum_filename):
        spectrum_filenames = [os.path.join(spectrum_filename, f) for f in os.listdir(spectrum_filename) if
                              f.endswith('.dat')]
        num_spectrum_files = len(spectrum_filenames)
        for file_index in range(num_spectrum_files):
            spectrum_filename_with_path = spectrum_filenames[file_index]
            i, e = spectrum_read(spectrum_filename_with_path)
            ivec.append(i)
            evec.append(e)
    elif os.path.isfile(spectrum_filename):
        num_spectrum_files = 1
        i, e = spectrum_read(spectrum_filename)
        ivec.append(i)
        evec.append(e)

    # Calculate and apply filter transmission
    for file_index in range(num_spectrum_files):
        mu_times_thickness = flat_filter(cfg.flat_filters, evec[file_index])
        ivec[file_index] = ivec[file_index] * np.exp(-mu_times_thickness)

    # Write spectrum
    if os.path.isdir(spectrum_filename):
        path_filtered = os.path.join(spectrum_filename, 'filtered')
    elif os.path.isfile(spectrum_filename):
        pathstr = os.path.dirname(spectrum_filename)
        path_filtered = os.path.join(pathstr, 'filtered')

    if not os.path.exists(path_filtered):
        os.makedirs(path_filtered)

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
