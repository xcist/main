# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import os
import numpy as np
from gecatsim.pyfiles.GetMu import GetMu

def catvoxel(configfilename, cfg, adjust, preadjust, silent):

    # SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    try:
        cfg.material_volumes
    except AttributeError:
        cfg.material_volumes = 0

    if cfg.material_volumes:
        print('Multiple volumes of material volume fractions will be produced.')
        if cfg.write_vp:
            print('Voxelized phantom files will be written.')
        else:
            print('Voxelized phantom files will not be written.')
    else:
        print(f"A single volume of attenuation coefficients will be produced.\nAttenuation coefficients will be determined at {cfg.make_img_kv} keV.")

    # PHANTOM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Materials, NumberOfMaterials = getattr(cfg, 'callback_getphantom')(cfg)

    # MATERIALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Mu = []
    for MaterialIndex in range(NumberOfMaterials):
        Mu.append(GetMu(Materials[MaterialIndex], cfg.make_img_kv))

    # Scale mu table in 1/mm and pass to C
    setattr(cfg, 'callback_setMaterial', lambda *args: None)  # Dummy function, not defined in provided code.
    getattr(cfg, 'callback_setMaterial')(NumberOfMaterials, 1, np.array(Mu)/10.0)

    try:
        cfg.Nx
    except AttributeError:
        cfg.Nx = cfg.recon_size
        cfg.Ny = cfg.recon_size
        cfg.Nz = cfg.recon_planes
        cfg.dx = cfg.recon_fov/cfg.recon_size
        cfg.dy = cfg.recon_fov/cfg.recon_size
        cfg.dz = cfg.recon_slice_thickness
        cfg.xoff = -cfg.recon_xcenter/cfg.dx + (cfg.Nx + 1) / 2
        cfg.yoff = -cfg.recon_ycenter/cfg.dy + (cfg.Ny + 1) / 2
        cfg.zoff = -cfg.recon_zcenter/cfg.dz + (cfg.Nz + 1) / 2

    print(f'Oversampling = {cfg.vol_os}')

    # VOLUMES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if cfg.material_volumes:
        # Make separate volumes of densities (in gm/cm^3) for each material
        print(f'Producing volume fraction volumes for {NumberOfMaterials} materials.')
    else:
        # Make just one volume (mu in cm^-1 at given kV)
        print('Producing a single volume of attenuation coefficients.')
        NumberOfMaterials = 1

    Volume = np.zeros((cfg.Nx, cfg.Ny, cfg.Nz, NumberOfMaterials))
    Volume = MakeAllVolumes(Volume, cfg, NumberOfMaterials, cfg.material_volumes)
    Volume = Volume.reshape(cfg.Nx, cfg.Ny, cfg.Nz, NumberOfMaterials)

    # WRITE FILEs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # Flip the volume top-to-bottom so it will be in the "file" format of decreasing Y value with increasing Y index.
    Volume_YFlipped = np.flip(Volume, axis=1)

    Directory, FileName, PhantomExtension = os.path.splitext(cfg.phantom_filename)
    if cfg.material_volumes and cfg.write_vp:
        # Write a .vp file and material density volume file(s)
        VoxelizedPhantomPathname = cfg.phantom_filename.replace(PhantomExtension, '.vp')
        print(f'Writing {VoxelizedPhantomPathname} and material volume fraction file(s)...')
        with open(VoxelizedPhantomPathname, 'w') as fid:
            fid.write(f'vp.n_materials = {NumberOfMaterials};\n')
            for MaterialIndex in range(NumberOfMaterials):
                VolumeExtension = f'VolumeFraction_{Materials[MaterialIndex]}'
                VolumePathname = VoxelizedPhantomPathname.replace('.vp', f'.{VolumeExtension}')
                Path, Filename, Extension = os.path.splitext(VolumePathname)
                VolumeFilename = Filename + Extension

                print(f'    writing material volume fraction file for {VolumeExtension}')
                with open(VolumePathname, 'wb') as f:
                    Volume_YFlipped[:, :, :, MaterialIndex].tofile(f)

                print(f'    writing data for {VolumeExtension} to {VoxelizedPhantomPathname}')
                fid.write(f'\nvp.mat_name{{{MaterialIndex + 1}}} = ''{Materials[MaterialIndex]}'';\n')
                fid.write(f'vp.volumefractionmap_filename{{{MaterialIndex + 1}}} = ''{VolumeFilename}'';\n')
                fid.write(f'vp.volumefractionmap_datatype{{{MaterialIndex + 1}}} = ''float'';\n')
                fid.write(f'vp.cols{{{MaterialIndex + 1}}} = {cfg.Nx};\n')
                fid.write(f'vp.rows{{{MaterialIndex + 1}}} = {cfg.Ny};\n')
                fid.write(f'vp.slices{{{MaterialIndex + 1}}} = {cfg.Nz};\n')
                fid.write(f'vp.x_size{{{MaterialIndex + 1}}} = {cfg.dx};\n')
                fid.write(f'vp.y_size{{{MaterialIndex + 1}}} = {cfg.dy};\n')
                fid.write(f'vp.z_size{{{MaterialIndex + 1}}} = {cfg.dz};\n')
                fid.write(f'vp.x_offset{{{MaterialIndex + 1}}} = {cfg.xoff};\n')
                fid.write(f'vp.y_offset{{{MaterialIndex + 1}}} = {cfg.yoff};\n')
                fid.write(f'vp.z_offset{{{MaterialIndex + 1}}} = {cfg.zoff};\n')

        print('... done writing files.')
    else:
        # Write attenuation coefficient volume to a .vol file
        VoxelizedPhantomPathname = cfg.phantom_filename.replace(PhantomExtension, '.mu')
        print(f'Writing attenuation coefficients volume to {VoxelizedPhantomPathname}')
        with open(VoxelizedPhantomPathname, 'wb') as f:
            Volume_YFlipped.tofile(f)

    # DONE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def MakeAllVolumes(Volume, cfg, NumberOfMaterials, MakeMaterialVolumes):
    switch_function = {
        'C_Projector_Analytic': 'Phantom_Analytic',
        'C_Projector_NCAT': 'Phantom_NCAT',
        'C_Projector_Polygon': 'Phantom_Polygonal'
    }

    FunctionName = switch_function.get(cfg.callback_projector, None)
    if FunctionName is None:
        raise ValueError(f"Error: CatVoxel not compatible with the projector_callback {cfg.callback_projector}")

    Volume = getattr(cfg, FunctionName)(
        Volume, cfg.Nx, cfg.xoff, cfg.dx,
        cfg.Ny, cfg.yoff, cfg.dy,
        cfg.Nz, cfg.zoff, cfg.dz,
        cfg.vol_os, NumberOfMaterials, MakeMaterialVolumes
    )

    return Volume