# Copyright 2025, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import os
from gecatsim.pyfiles.catvoxel import catvoxel
import numpy as np

def Phantom_Analytic_to_Voxelized_to_VolumesOfMassConcentrations(cfg):

    directory, filename = os.path.split(cfg.phantom_filename)
    name, extension = os.path.splitext(filename)
    if extension not in ['.pp', '.ppm']:
        raise ValueError(f'Unknown phantom file format: {cfg.phantom_filename}')

    voxelized_phantom_filename = cfg.phantom_filename.replace(extension, '.vp')

    if not os.path.exists(voxelized_phantom_filename) or \
       os.path.getmtime(cfg.phantom_filename) > os.path.getmtime(voxelized_phantom_filename) or \
       (hasattr(cfg, 'force_phantom_conversion') and cfg.force_phantom_conversion):

        print('Voxelizing analytic phantom...')

        print('Setting phantom matrix per cfg.phantom_samples_*')
        adjust = [
            'cfg.material_volumes = 1;',
            f'cfg.Nx = {cfg.phantom_samples_xy};',
            f'cfg.dx = {cfg.phantom_samples_voxelsize};',
            f'cfg.xoff = {(cfg.phantom_samples_xy + 1) / 2};',
            f'cfg.Ny = {cfg.phantom_samples_xy};',
            f'cfg.dy = {cfg.phantom_samples_voxelsize};',
            f'yoff = {(cfg.phantom_samples_xy + 1) / 2};',
            f'cfg.Nz = {cfg.phantom_samples_z};',
            f'cfg.dz = {cfg.phantom_samples_voxelsize};',
            f'zoff = {(cfg.phantom_samples_z + 1) / 2};'
        ]

        cfg.write_vp = 1
        Volume, MaterialList = catvoxel([], cfg, adjust)

        print('... done voxelizing analytic phantom.')

    cfg.phantom_filename = voxelized_phantom_filename
    # Phantom_Voxelized_to_VolumesOfMassConcentrations() # restructured to phantom_Voxelized
