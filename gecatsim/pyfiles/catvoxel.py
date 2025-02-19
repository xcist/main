# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license
"""
-----------------------------------------------------------------------
GE proprietary and confidential
CatSim Detector Model for Revolution 

Authors
  Mingye Wu (GE Global Research)

Aim
  Returns the detector geometry parameters.
  This function itself calculates the module coordinates and UV vectors.
  It calls Pack_Module.m to obtain the parameters of cells and
  oversamples, and calls Ideal_Grid.m to obtain the grid data.

Inputs
  cfg.sid : Source-to-isocenter distance (mm)
  cfg.sdd : Source-to-detector distance (mm)
  cfg.callback_detector : Build a curved fs-centered segmented detector for SVCT
  cfg.col_size : Detector column size or x-pitch (mm)
  cfg.row_size : Detector row size or z-pitch (mm)
  cfg.col_cast : Cast reflector thickness at column direction, in mm
  cfg.row_cast : Cast reflector thickness at raw direction, in mm
  cfg.mods_per_det_x : Number of modules per detector
  cfg.mods_per_det_z : Number of modules per detector
  cfg.mod_gap_u : Inter-module gap at column direction (in mm)
  cfg.mod_gap_v : Inter-module gap at row direction (in mm)
  cfg.col_offset : In-plane detector offset (in pixel)
  cfg.mod_row_offset : [n_modules] vector with longitudinal detector offset for each module (in pixel)
  cfg.cols_per_mod : Number of detector columns per module
  cfg.rows_per_mod : Number of detector rows per module
  cfg.plate_thickness : 1D post-collimator plate thickness (in mm)
  cfg.plate_height : 1D post-collimator plate height (in mm)
  cfg.plate_airgap : Distance from bottom of post-collimator plate to scintillator surface (in mm)
  cfg.detector_depth : Depth of the conversion layer (mm)
  cfg.col_oversample : Detector column normal oversampling. Need to be even number for collimator plate at cell center.
  cfg.row_oversample : Detector row oversampling
  cfg.col_crosstalk : fractional crosstalk to each neighboring column
                      eg: 0.02 results in 96% weight for current column
                      and 2% for each neighboring column
  cfg.row_crosstalk : fractional crosstalk to each neighboring row
                      eg: 0.02 results in 96% weight for current row
                      and 2% for each neighboring row

Outputs
  det : structure containing the following fields :
    n_modules     : number of detector modules
    modtypes      : vector of size n_modules with the types of each module :
                    this index is one-based
    modcoords     : [3 x n_modules] array with the module center
                    xyz coordinates (in mm)
    uvecs         : [3 x n_modules] array with the transversal unit vectors
    vvecs         : [3 x n_modules] array with the longitudinal unit vectors
    total_n_cells : total number of cells in the detector
    n_moddefs     : number of different module types
    n_cells       : [n_moddefs] vector with number of cells in each module type
    cellcoords    : [2 x max(n_cells) x n_moddefs] array of (u,v) coords for each module type
                    cellcoords are relative to module center, uvec and vvec, in mm
    n_samples     : [n_moddefs] vector with number of samples per cell
    samplecoords  : [2 x max(n_samples) x n_moddefs] array of (u,v,w)
                    u and v are the sample coords relative to cell center
                    uvec and vvec, in mm. The sample locations may include x-ray cross-talk samples in neighboring cells.
    weights       : [max(n_samples) x n_moddefs] array of relative sample weights. The weights have to sum up to one.
    activearea    : [n_moddefs] vector with activearea of the cells
    width         : [n_moddefs] vector with module widths, serves as bounding box for all cell samples
    height        : [n_moddefs] vector with module heights, serves as bounding box for all cell samples
    startindices  : [n_modules] array containing 0-based indices into the view for the first pixel in each module

Note:
   Support only one module type
   Anti-Scatter-Grid is implemented in catsim.m.
-----------------------------------------------------------------------
"""
import os
import numpy as np
from gecatsim.pyfiles.GetMu import GetMu
from gecatsim.pyfiles.CommonTools import *

def catvoxel(configfilename, cfg, adjust=None, preadjust=None, silent=False):

    cfg.phantom.filename = my_path.find('phantom', cfg.phantom.filename, '')
    filepath, vpFilename = os.path.split(cfg.phantom.filename)
    with open(cfg.phantom.filename) as fin:
        vp = json.load(fin)

    ###----------- pass material Mu to C
    Materials = vp['mat_name']
    NumberOfMaterials = vp['n_materials']

    # MATERIALS
    Mu = []
    for MaterialIndex in range(NumberOfMaterials):
        Mu.append(GetMu(Materials[MaterialIndex], cfg.make_img_kv))

    # Scale mu table in 1/mm and pass to C
    cfg.phantom.setMaterial(NumberOfMaterials, 1, np.array(Mu) / 10.0)

    try:
        cfg.Nx
    except AttributeError:
        cfg.Nx = cfg.recon_size
        cfg.Ny = cfg.recon_size
        cfg.Nz = cfg.recon_planes
        cfg.dx = cfg.recon_fov / cfg.recon_size
        cfg.dy = cfg.recon_fov / cfg.recon_size
        cfg.dz = cfg.recon_slice_thickness
        cfg.xoff = -cfg.recon_xcenter / cfg.dx + (cfg.Nx + 1) / 2
        cfg.yoff = -cfg.recon_ycenter / cfg.dy + (cfg.Ny + 1) / 2
        cfg.zoff = -cfg.recon_zcenter / cfg.dz + (cfg.Nz + 1) / 2
        print('Phantom setup updated:')
        PrintPhantomSetup(cfg)

    print(f'Oversampling = {cfg.vol_os}')

    # VOLUMES
    if cfg.material_volumes:
        print(f'Producing volume fraction volumes for {NumberOfMaterials} materials.')
    else:
        print('Producing a single volume of attenuation coefficients.')
        NumberOfMaterials = 1

    Volume = np.zeros((cfg.Nx, cfg.Ny, cfg.Nz, NumberOfMaterials))
    Volume = MakeAllVolumes(Volume, cfg, NumberOfMaterials, cfg.material_volumes)

    # WRITE FILEs

    # Flip the volume top-to-bottom so it will be in the "file" format of decreasing Y value with increasing Y index.
    Volume_YFlipped = np.flip(Volume, axis=1)

    Directory, PhantomExtension = os.path.splitext(cfg.phantom.filename)
    FileName = os.path.basename(Directory)
    if cfg.material_volumes and cfg.write_vp:
        # Write a .vp file and material density volume file(s)
        VoxelizedPhantomPathname = cfg.phantom.filename.replace(PhantomExtension, '.vp')
        print(f'Writing {VoxelizedPhantomPathname} and material volume fraction file(s)...')
        with open(VoxelizedPhantomPathname, 'w') as fid:
            fid.write(f'vp.n_materials = {NumberOfMaterials};\n')
            for MaterialIndex in range(NumberOfMaterials):
                VolumeExtension = f'VolumeFraction_{Materials[MaterialIndex]}'
                VolumePathname = VoxelizedPhantomPathname.replace('.vp', f'.{VolumeExtension}')
                Path, Filename = os.path.splitext(VolumePathname)
                VolumeFilename = Filename + PhantomExtension

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
        VoxelizedPhantomPathname = cfg.phantom.filename.replace(PhantomExtension, '.mu')
        print(f'Writing attenuation coefficients volume to {VoxelizedPhantomPathname}')
        with open(VoxelizedPhantomPathname, 'wb') as f:
            Volume_YFlipped.tofile(f)

    # DONE


def MakeAllVolumes(Volume, cfg, NumberOfMaterials, MakeMaterialVolumes):
    switch_function = {
        'C_Projector_Analytic': 'Phantom_Analytic',
        'C_Projector_NCAT': 'Phantom_NCAT',
        'C_Projector_Polygon': 'Phantom_Polygonal'
    }

    FunctionName = switch_function.get(cfg.phantom.projectorCallback, None)
    if FunctionName is None:
        raise ValueError(f"Error: CatVoxel not compatible with the projector_callback {cfg.phantom.projectorCallback}")

    Volume = getattr(cfg.phantom, FunctionName)(
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