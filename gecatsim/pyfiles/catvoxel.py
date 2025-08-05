# Copyright 2025, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license
"""
-----------------------------------------------------------------------
convert analytic, ncat, or polygonal phantoms to voxelized phantom

Authors
  Mingye Wu (GE Global Research)
  Jiayong Zhang

-----------------------------------------------------------------------
"""

import os
import numpy as np
import gecatsim as xc
import json
from collections import defaultdict
from gecatsim.pyfiles.GetMu import GetMu
from gecatsim.pyfiles.CommonTools import *

def dumpjson(filename, jsondict):
    out_file = open(filename, "w") 
    json.dump(jsondict, out_file, indent = 4) 
    out_file.close()

def catvoxel(cfg):

    if not hasattr(cfg, "material_volumes"):
        cfg.material_volumes = 0
    if cfg.material_volumes:
        print('Multiple volumes of material volume fractions will be produced.')
        if cfg.write_vp:
            print('Voxelized phantom files will be written.')
        else:
            print('Voxelized phantom files will not be written.')
    else:
        if hasattr(cfg, "make_img_kv"):
            cfg.spec.Evec = np.array(cfg.make_img_kv)
        else:
            cfg.spec.Evec = np.array([120])

        print(f"A single volume of attenuation coefficients will be produced.\nAttenuation coefficients will be determined at {cfg.make_img_kv} keV.")

    # PHANTOM 
    cfg = feval(cfg.phantom.callback, cfg)
    Materials = cfg.phantom.Materials
    NumberOfMaterials = cfg.phantom.numberOfMaterials

    cfg.Nx = cfg.recon.imageSize
    cfg.Ny = cfg.recon.imageSize
    cfg.Nz = cfg.recon.sliceCount
    cfg.dx = cfg.recon.fov / cfg.recon.imageSize
    cfg.dy = cfg.recon.fov / cfg.recon.imageSize
    cfg.dz = cfg.recon.sliceThickness
    cfg.xoff = -cfg.recon.centerOffset[0] / cfg.dx + (cfg.Nx + 1) / 2
    cfg.yoff = -cfg.recon.centerOffset[1] / cfg.dy + (cfg.Ny + 1) / 2
    cfg.zoff = -cfg.recon.centerOffset[2] / cfg.dz + (cfg.Nz + 1) / 2
    print('Phantom setup updated:')
    PrintPhantomSetup(cfg)

    print(f'Oversampling = {cfg.vol_os}')

    # VOLUMES
    if cfg.material_volumes:
        print(f'Producing volume fraction volumes for {NumberOfMaterials} materials.')
    else:
        print('Producing a single volume of attenuation coefficients.')
        NumberOfMaterials = 1

    Volume = np.zeros((NumberOfMaterials, cfg.Nz, cfg.Ny, cfg.Nx), dtype=np.single)
    Volume = MakeAllVolumes(Volume, cfg, NumberOfMaterials, cfg.material_volumes)
    Volume = np.transpose(Volume, axes=[3, 2, 1, 0])

    # Flip the volume top-to-bottom so it will be in the "file" format of decreasing Y value with increasing Y index.
    Volume_YFlipped = np.flip(Volume, axis=1)

    Directory, PhantomExtension = os.path.splitext(cfg.phantom.filename)
    FileName = os.path.basename(Directory)

    # to avoid the case when phantom is in python library, the output will be saved there
    base_phantom_name = os.path.basename(cfg.phantom.filename)
    if cfg.material_volumes and cfg.write_vp:
        vp_json = defaultdict(list)
        nonzero_materials = 0
        # Write a .vp file and material density volume file(s)
        VoxelizedPhantomPathname = base_phantom_name.replace(PhantomExtension, '.json')
        print(f'Writing {VoxelizedPhantomPathname} and material volume fraction file(s)...')
        # this line needs to be here so that in json file, the position of n_materials is right
        vp_json['n_materials'] = 0
        for MaterialIndex in range(NumberOfMaterials):

            # ignore this material if all zero
            if np.sum(Volume[:,:,:,MaterialIndex]) < 1.E-4: continue

            nonzero_materials += 1

            VolumeExtension = f'VolumeFraction_{Materials[MaterialIndex]}'
            VolumePathname = VoxelizedPhantomPathname.replace('.json', f'.{VolumeExtension}')
            Path, Filename = os.path.splitext(VolumePathname)
            VolumeFilename = Filename + PhantomExtension
            print(f'    writing material volume fraction file for {VolumeExtension}')

            tmpArr = np.transpose(Volume_YFlipped[:, :, :, MaterialIndex], axes=[2, 1, 0])
            xc.rawwrite(VolumePathname, np.copy(tmpArr, order='C'))

            print(f'    writing data for {VolumeExtension} to {VoxelizedPhantomPathname}')
            vp_json['mat_name'].append(f'{Materials[MaterialIndex]}')
            vp_json['volumefractionmap_filename'].append(VolumePathname)
            vp_json['volumefractionmap_datatype'].append('float')
            vp_json['cols'].append(cfg.Nx)
            vp_json['rows'].append(cfg.Ny)
            vp_json['slices'].append(cfg.Nz)
            vp_json['x_size'].append(cfg.dx)
            vp_json['y_size'].append(cfg.dy)
            vp_json['z_size'].append(cfg.dz)
            vp_json['x_offset'].append(cfg.xoff)
            vp_json['y_offset'].append(cfg.yoff)
            vp_json['z_offset'].append(cfg.zoff)

        vp_json['n_materials'] = nonzero_materials
        dumpjson(VoxelizedPhantomPathname, vp_json)
        print('... done writing files.')
    else:
        # Write attenuation coefficient volume to a .vol file
        VoxelizedPhantomPathname = base_phantom_name.replace(PhantomExtension, '.mu')
        print(f'Writing attenuation coefficients volume to {VoxelizedPhantomPathname}')
        with open(VoxelizedPhantomPathname, 'wb') as f:
            Volume_YFlipped.tofile(f)


def MakeAllVolumes(Volume, cfg, NumberOfMaterials, MakeMaterialVolumes):
    switch_function = {
        'C_Projector_Analytic': 'C_Volume_Get',
        'C_Projector_NCAT': 'C_Volume_NCAT_Get',
        'C_Projector_Polygon': 'C_Volume_Poly_Get'
    }

    FunctionName = switch_function.get(cfg.phantom.projectorCallback, None)
    if FunctionName is None:
        raise ValueError(f"Error: CatVoxel not compatible with the projector_callback {cfg.phantom.projectorCallback}")

    Volume = feval(FunctionName, cfg, 
        Volume, cfg.Nx, cfg.xoff, cfg.dx,
        cfg.Ny, cfg.yoff, cfg.dy,
        cfg.Nz, cfg.zoff, cfg.dz,
        cfg.vol_os, NumberOfMaterials, MakeMaterialVolumes
    )

    return Volume

def PrintPhantomSetup(cfg):
    print(f'Number of X voxels: cfg.Nx   = {cfg.Nx}')
    print(f'Number of Y voxels: cfg.Ny   = {cfg.Ny}')
    print(f'Number of Z voxels: cfg.Nz   = {cfg.Nz}')
    print(f'X voxel size: cfg.dx  = {cfg.dx}')
    print(f'Y voxel size: cfg.dy  = {cfg.dy}')
    print(f'Z voxel size: cfg.dz  = {cfg.dz}')
    print(f'X voxel offset: cfg.xoff = {cfg.xoff}')
    print(f'Y voxel offset: cfg.yoff = {cfg.yoff}')
    print(f'Z voxel offset: cfg.zoff = {cfg.zoff}')

