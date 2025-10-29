# Copyright 2025, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

# -----------------------------------------------------------------------
# Aim:
#     This script runs `catvoxel` if needed to produce a set of material volumes,
#     and invokes `dose_map_v`.
#
# Input:
#     Assumes:
#         - A valid `cfg` structure.
#         - `cfg.phantom_filename` contains a valid path to an analytic phantom.
#
# Output:
#     - A voxelized phantom file equivalent to the analytic phantom described in
#       `cfg.phantom_filename` at entry.
#     - `cfg.phantom_filename` is updated with the name of the voxelized phantom.
# -----------------------------------------------------------------------

import os
from gecatsim.pyfiles.catvoxel import catvoxel

def Phantom_Analytic_to_Voxelized_to_VolumesOfMassConcentrations(cfg):
    phantom_file = cfg.phantom.filename
    directory, filename = os.path.split(phantom_file)
    name, extension = os.path.splitext(filename)

    if extension not in ['.pp', '.ppm']:
        raise ValueError(f'Unknown phantom file format: {phantom_file}')

    voxelized_phantom_filename = phantom_file.replace(extension, '.vp')

    if not os.path.exists(voxelized_phantom_filename) or \
       os.path.getmtime(phantom_file) > os.path.getmtime(voxelized_phantom_filename) or \
       (hasattr(cfg, 'force_phantom_conversion') and cfg.force_phantom_conversion):

        print('Voxelizing analytic phantom...')
        print('Setting phantom matrix per cfg.phantom_samples_*')

        cfg.material_volumes = 1
        cfg.Nx = cfg.phantom.samples_xy
        cfg.dx = cfg.phantom.samples_voxelsize
        cfg.xoff = (cfg.Nx + 1) / 2
        cfg.Ny = cfg.phantom.samples_xy
        cfg.dy = cfg.phantom.samples_voxelsize
        cfg.yoff = (cfg.Ny + 1) / 2
        cfg.Nz = cfg.phantom.samples_z
        cfg.dz = cfg.phantom.samples_voxelsize
        cfg.zoff = (cfg.Nz + 1) / 2

        cfg.write_vp = 1
        if not hasattr(cfg, "make_img_kv"):
            cfg.make_img_kv = 120  # Default value

        Volume, MaterialList = catvoxel(cfg)

        print('... done voxelizing analytic phantom.')

    cfg.phantom_filename = voxelized_phantom_filename
