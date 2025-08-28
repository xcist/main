import os
from gecatsim.pyfiles.catvoxel import catvoxel
import numpy as np

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
        print('Setting phantom matrix per cfg.phantom.samples_*')

        adjust = [
            'cfg.material_volumes = 1;',
            f'cfg.Nx = {cfg.phantom.samples_xy};',
            f'cfg.dx = {cfg.phantom.samples_voxelsize};',
            f'cfg.xoff = {(cfg.phantom.samples_xy + 1) / 2};',
            f'cfg.Ny = {cfg.phantom.samples_xy};',
            f'cfg.dy = {cfg.phantom.samples_voxelsize};',
            f'yoff = {(cfg.phantom.samples_xy + 1) / 2};',
            f'cfg.Nz = {cfg.phantom.samples_z};',
            f'cfg.dz = {cfg.phantom.samples_voxelsize};',
            f'zoff = {(cfg.phantom.samples_z + 1) / 2};'
        ]

        cfg.write_vp = 1
        cfg.material_volumes = 1  # Ensure this is set
        if not hasattr(cfg, "make_img_kv"):
            cfg.make_img_kv = 120  # Default value to avoid AttributeError

        Volume, MaterialList = catvoxel(cfg)

        print('... done voxelizing analytic phantom.')

    cfg.phantom.filename = voxelized_phantom_filename
