# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import os
import numpy as np
from gecatsim.pyfiles.GetMu import GetMu
import gecatsim.pyfiles.CatSim as xc

def catvoxel(configfilename, cfg, adjust, preadjust, silent):

    try:
        cfg['material_volumes']
    except KeyError:
        cfg['material_volumes'] = 0

    # if cfg['material_volumes']:
    #     Headline(['Multiple volumes of material volume fractions will be produced.'], Verbose.Phantom)
    #     if cfg['write_vp']:
    #         Headline(['Voxelized phantom files will be written.'], Verbose.Phantom)
    #     else:
    #         Headline(['Voxelized phantom files will not be written.'], Verbose.Phantom)
    # else:
    #     Headline(['A single volume of attenuation coefficients will be produced.',
    #               f"Attenuation coefficients will be determined at {cfg['make_img_kv']} keV."], Verbose.Phantom)

    # PHANTOM
    Materials, NumberOfMaterials = cfg.callback_getphantom

    # MATERIALS
    Mu = []
    for MaterialIndex in range(NumberOfMaterials):
        Mu.append(GetMu(Materials[MaterialIndex], cfg['make_img_kv']))

    # Scale mu table in 1/mm and pass to C
    mu = np.array(Mu) / 10.0
    cfg.callback_setMaterial(NumberOfMaterials, 1, mu)

    try:
        cfg['Nx']
    except KeyError:
        cfg['Nx'] = cfg['recon_size']
        cfg['Ny'] = cfg['recon_size']
        cfg['Nz'] = cfg['recon_planes']
        cfg['dx'] = cfg['recon_fov'] / cfg['recon_size']
        cfg['dy'] = cfg['recon_fov'] / cfg['recon_size']
        cfg['dz'] = cfg['recon_slice_thickness']
        cfg['xoff'] = -cfg['recon_xcenter'] / cfg['dx'] + (cfg['Nx'] + 1) / 2
        cfg['yoff'] = -cfg['recon_ycenter'] / cfg['dy'] + (cfg['Ny'] + 1) / 2
        cfg['zoff'] = -cfg['recon_zcenter'] / cfg['dz'] + (cfg['Nz'] + 1) / 2
        # if Verbose.Phantom:
        #     PrintPhantomSetup(cfg)

    # Headline([f'Oversampling = {cfg["vol_os"]}'], Verbose.Phantom)

    # VOLUMES
    if not cfg['material_volumes']:
        # Headline([f'Producing volume fraction volumes for {NumberOfMaterials} materials.'], Verbose.Phantom)
    # else:
        # Headline(['Producing a single volume of attenuation coefficients.'], Verbose.Phantom)
        NumberOfMaterials = 1

    Volume = np.zeros((cfg['Nx'], cfg['Ny'], cfg['Nz'], NumberOfMaterials))
    Volume = MakeAllVolumes(Volume, cfg, NumberOfMaterials, cfg['material_volumes'])
    Volume = Volume.reshape((cfg['Nx'], cfg['Ny'], cfg['Nz'], NumberOfMaterials))

    # WRITE FILES
    Volume_YFlipped = np.flip(Volume, 1)

    Directory, FileName, PhantomExtension = os.path.splitext(cfg['phantom_filename'])
    if cfg['material_volumes']:
        if cfg['write_vp']:
            VoxelizedPhantomPathname = cfg['phantom_filename'].replace(PhantomExtension, '.vp')
            # Headline([f'Writing {VoxelizedPhantomPathname} and material volume fraction file(s)...'], Verbose.Phantom)
            with open(VoxelizedPhantomPathname, 'w') as fid:
                fid.write(f'vp.n_materials = {NumberOfMaterials};\n')
                for MaterialIndex in range(NumberOfMaterials):
                    VolumeExtension = f'VolumeFraction_{Materials[MaterialIndex]}'
                    VolumePathname = VoxelizedPhantomPathname.replace('.vp', f'.{VolumeExtension}')
                    Path, Filename, Extension = os.path.splitext(VolumePathname)
                    VolumeFilename = f'{Filename}{Extension}'
                    # Headline([f'    writing material volume fraction file for {VolumeExtension}'], Verbose.Phantom)
                    xc.rawwrite(VolumePathname, Volume_YFlipped[:, :, :, MaterialIndex].astype(np.float32))
                    # Headline([f'    writing data for {VolumeExtension} to {VoxelizedPhantomPathname}'], Verbose.Phantom)
                    fid.write(f"\nvp.mat_name{{{MaterialIndex + 1}}} = '{Materials[MaterialIndex]}';\n")
                    fid.write(f"vp.volumefractionmap_filename{{{MaterialIndex + 1}}} = '{VolumeFilename}';\n")
                    fid.write(f"vp.volumefractionmap_datatype{{{MaterialIndex + 1}}} = 'float';\n")
                    fid.write(f"vp.cols{{{MaterialIndex + 1}}} = {cfg['Nx']};\n")
                    fid.write(f"vp.rows{{{MaterialIndex + 1}}} = {cfg['Ny']};\n")
                    fid.write(f"vp.slices{{{MaterialIndex + 1}}} = {cfg['Nz']};\n")
                    fid.write(f"vp.x_size{{{MaterialIndex + 1}}} = {cfg['dx']};\n")
                    fid.write(f"vp.y_size{{{MaterialIndex + 1}}} = {cfg['dy']};\n")
                    fid.write(f"vp.z_size{{{MaterialIndex + 1}}} = {cfg['dz']};\n")
                    fid.write(f"vp.x_offset{{{MaterialIndex + 1}}} = {cfg['xoff']};\n")
                    fid.write(f"vp.y_offset{{{MaterialIndex + 1}}} = {cfg['yoff']};\n")
                    fid.write(f"vp.z_offset{{{MaterialIndex + 1}}} = {cfg['zoff']};\n")

    else:
        VoxelizedPhantomPathname = cfg['phantom_filename'].replace(PhantomExtension, '.mu')
        # Headline([f'Writing attenuation coefficients volume to {VoxelizedPhantomPathname}'], Verbose.Phantom)
        xc.rawwrite(VoxelizedPhantomPathname, Volume_YFlipped.astype(np.float32))

    # DONE

def MakeAllVolumes(Volume, cfg, NumberOfMaterials, MakeMaterialVolumes):
    if cfg['callback_projector'] == 'C_Projector_Analytic':
        FunctionName = 'C_Volume_Get'
    elif cfg['callback_projector'] == 'C_Projector_NCAT':
        FunctionName = 'C_Volume_NCAT_Get'
    elif cfg['callback_projector'] == 'C_Projector_Polygon':
        FunctionName = 'C_Volume_NCAT_Get'
    else:
        raise ValueError(f"Error: CatVoxel not compatible with the projector_callback {cfg['callback_projector']}")

    Volume = eval(FunctionName)(Volume, cfg['Nx'], cfg['xoff'], cfg['dx'],
                                cfg['Ny'], cfg['yoff'], cfg['dy'],
                                cfg['Nz'], cfg['zoff'], cfg['dz'],
                                cfg['vol_os'], NumberOfMaterials, MakeMaterialVolumes)
    return Volume

def PrintPhantomSetup(cfg):
    pass
    # global Verbose
    #
    # if Verbose.Phantom:
    #     Headline([' ',
    #               f'Number of X voxels: cfg.Nx   = {cfg["Nx"]}',
    #               f'Number of Y voxels: cfg.Ny   = {cfg["Ny"]}',
    #               f'Number of Z voxels: cfg.Nz   = {cfg["Nz"]}',
    #               ' ',
    #               f'      X voxel size: cfg.dNx  = {cfg["dx"]}',
    #               f'      Y voxel size: cfg.dNy  = {cfg["dy"]}',
    #               f'      Z voxel size: cfg.dNz  = {cfg["dz"]}',
    #               ' ',
    #               f'    X voxel offset: cfg.xoff = {cfg["xoff"]}',
    #               f'    Y voxel offset: cfg.yoff = {cfg["yoff"]}',
    #               f'    Z voxel offset: cfg.zoff = {cfg["zoff"]}',
    #               ' '], Verbose.Phantom)