# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import os
from gecatsim.pyfiles.CommonTools import *

def phantom_dynamic_update(cfg, iscatscatter=False):
    # if not iscatscatter:
    #     iscatscatter = False
    #
    # if cfg.dynamic_phantom:
    #     need_update = False
    #     if cfg.current_phase == 0:  # need to initialize
    #         cfg.phase_start_time = cfg.start_time
    #         cfg.current_phase = 1
    #
    #         # If we are simulating only a subset of the views, we need to change the phase until it catches up to
    #         # the first view we are simulating
    #         while cfg.time - cfg.phase_start_time >= cfg.phase_length(cfg.current_phase):
    #             cfg.phase_start_time += cfg.phase_length(cfg.current_phase)
    #             cfg.current_phase = (cfg.current_phase % cfg.num_phases) + 1
    #             need_update = True
    #
    #     while cfg.time - cfg.phase_start_time >= cfg.phase_length(cfg.current_phase):
    #         cfg.phase_start_time += cfg.phase_length(cfg.current_phase)
    #         cfg.current_phase = (cfg.current_phase % cfg.num_phases) + 1
    #         need_update = True
    #
    #     cfg.phantom_update_next_subview = 0
    #     time_next_subview = cfg.time + cfg.subviewtime
    #     if cfg.duty_ratio < 1 and cfg.subviewindex == cfg.subviews_per_view:
    #         time_next_subview += (1 - cfg.duty_ratio) * cfg.ViewTime
    #     if time_next_subview - cfg.phase_start_time >= cfg.phase_length(cfg.current_phase):
    #         cfg.phantom_update_next_subview = 1
    #
    #     if need_update:
    #         cfg.phantom_filename = cfg.phantom_filename_template.replace('&', str(cfg.current_phase))
    #         if iscatscatter:
    #             if 'vp' not in cfg.phantom_filename:
    #                 phantom_analytic_to_voxelized_to_volumes_of_mass_concentrations()
    #             else:
    #                 phantom_voxelized_to_volumes_of_mass_concentrations()
    #         else:
    #             materials, n_materials = cfg.callback_getphantom

    pass
