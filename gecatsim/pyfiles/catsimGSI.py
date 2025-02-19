# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

"""
Aim
    Simulate a realistic CT scan and write out projection data with GSI mode.
    It's also a wrapper of catsim.m.

For info on input and output arguments for the CatSim family applications, type
    help Cat_Common_IO
"""
import os
import numpy as np
import gecatsim as xc
from gecatsim.pyfiles.CommonTools import rawread
from gecatsim.pyfiles.PrepView import prep_view

def catsimGSI(cfg, adjust, preadjust, silent):
    kVp_0 = cfg.kVp  # note the initial kVp
    ViewTime = cfg.rotation_period / cfg.views_per_rotation
    start_time = cfg.start_time + (cfg.start_view - 1.5) * ViewTime

    # GSI parameters
    n_subphase = len(cfg.duty_cycle_per_view)
    n_phase = len(cfg.subphase_grouping)
    cfg.duty_cycle_per_view = np.array(cfg.duty_cycle_per_view) / sum(cfg.duty_cycle_per_view)  # Convert to NumPy array
    duty_cycle = cfg.duty_cycle_per_view

    airscan = [None] * n_phase
    offsetscan = [None] * n_phase
    view_projection = [None] * n_phase

    # process catsim simulation for each sub-phase
    print("\nGenerating projections: Beginning sub-phase loop...\n")
    for subphase in range(n_subphase):
        tmp_cfg = cfg

        # file basename
        tmp_cfg.results_basename = f"{cfg.results_basename}_subphase_{subphase + 1}"

        # kVp (spectrum) and mA
        kVp = cfg.kVp_cycle_per_view[subphase]
        tmp_cfg.number_Ebins = int(cfg.number_Ebins * kVp / kVp_0)
        tmp_cfg.spectrum_filename = cfg.spectrum_filename.replace(f"_{kVp_0}_", f"_{kVp}_")
        tmp_cfg.mA = cfg.mA_cycle_per_view[subphase]

        # duty ratio
        tmp_cfg.duty_ratio = duty_cycle[subphase]

        # start time adjust
        if tmp_cfg.start_time_GSI_mode == 1:
            tmp_cfg.start_time = start_time + (sum(duty_cycle[:subphase]) + 0.5 * duty_cycle[subphase]) * ViewTime
        else:
            tmp_cfg.start_time = start_time + 0.5 * ViewTime

        # other adjusts
        tmp_cfg.convert_to_prep = 0

        # process Duty Cycle Simulation
        print(
            f"\nCalling catsim for sub-phase {subphase + 1} of {n_subphase}, at {kVp:.0f} kVp, {tmp_cfg.mA:.0f} mA, {100 * tmp_cfg.duty_ratio:.0f} % duty cycle.\n")
        xc.CatSim([], tmp_cfg, [], [], 1)
        print(f"\nProjection generation for sub-phase {subphase + 1} of {n_subphase} complete.\n")

    print("\nProjection generation complete.\n")

    # read and process scan data by Duty Cycle Simulations
    if not cfg.this_is_an_airscan and not cfg.this_is_an_offsetscan:
        print("\nProcessing projections: Beginning phase loop...\n")
        for phase in range(n_phase):
            print(f"\nProcessing phase {phase + 1} of {n_phase}.\n")
            # Initialize arrays to accumulate data
            airscan[phase] = np.zeros((cfg.total_n_cells, 1, cfg.airscan_total_n_views))
            offsetscan[phase] = np.zeros((cfg.total_n_cells, 1, cfg.offsetscan_total_n_views))
            view_projection[phase] = np.zeros((cfg.total_n_cells, 1, cfg.total_n_views))

            # read scan data
            for subphase in cfg.subphase_grouping[phase]:
                print(f"\nReading data for sub-phase {subphase} of {n_subphase}.\n")
                try:
                    airscan_tmp = rawread(f"{cfg.results_basename}_subphase_{subphase}.air",
                                          (cfg.total_n_cells, 1, cfg.airscan_total_n_views), 'float')
                    offsetscan_tmp = rawread(f"{cfg.results_basename}_subphase_{subphase}.offset",
                                             (cfg.total_n_cells, 1, cfg.offsetscan_total_n_views), 'float')
                    view_projection_tmp = rawread(f"{cfg.results_basename}_subphase_{subphase}.scan",
                                                  (cfg.total_n_cells, 1, cfg.total_n_views), 'float')

                    airscan[phase] += airscan_tmp
                    offsetscan[phase] += offsetscan_tmp
                    view_projection[phase] += view_projection_tmp
                except FileNotFoundError as e:
                    print(f"Warning: {e}. Skipping this sub-phase.")

            # generate prep
            extension = '.scan'
            if cfg.convert_to_prep_per_phase == 1:
                print("\nPrepping views.\n")
                for kk in range(cfg.start_view, cfg.stop_view + 1):
                    view_projection[phase][:, :, kk] = prep_view(cfg)
                extension = '.prep'

            # write to files
            print("\nWriting files.\n")
            with open(f"{cfg.results_basename}_phase_{phase + 1}.air", 'wb') as fid:
                fid.write(airscan[phase].astype('float32').tobytes())

            with open(f"{cfg.results_basename}_phase_{phase + 1}.offset", 'wb') as fid:
                fid.write(offsetscan[phase].astype('float32').tobytes())

            with open(f"{cfg.results_basename}_phase_{phase + 1}{extension}", 'wb') as fid:
                fid.write(view_projection[phase].astype('float32').tobytes())

            print(f"\nProjection processing for phase {phase + 1} of {n_phase} complete.\n")
        print("\nProjection processing complete.\n")