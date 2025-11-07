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
import numpy as np
from gecatsim.pyfiles.CommonTools import *
from gecatsim.pyfiles.Detector_Pack_Module import detector_pack_module

def detector_thirdgen_curved_multi_z_modules(cfg):

    print("\nComputing a THIRD-GENERATION CURVED detector with MULTIPLE MODULES IN Z...\n")

    # SHORTCUTS
    sid = cfg.scanner.sid
    sdd = cfg.scanner.sdd

    colsize = cfg.scanner.detectorColSize
    rowsize = cfg.scanner.detectorRowSize
    n_rows = cfg.scanner.detectorRowsPerMod
    n_cols = cfg.scanner.detectorColsPerMod

    n_modules_x = cfg.scanner.modsPerDetX
    n_modules_z = cfg.scanner.modsPerDetZ
    n_modules = int(n_modules_x * n_modules_z)

    mod_gap_u = cfg.scanner.modGapU
    mod_gap_v = cfg.scanner.modGapV

    mod_offset_w = cfg.scanner.modOffsetW
    mod_rotation_u = cfg.scanner.modRotationU

    # Module alphas in X-Y plan
    modwidth_x = n_cols * colsize + cfg.scanner.colPackGap * (cfg.scanner.colPacksPerMod - 1)
    dalpha_mod = 2 * np.arctan(modwidth_x * 0.5 / sdd)
    edge_sdd = np.sqrt((modwidth_x * 0.5) ** 2 + sdd ** 2)
    dalpha_modgap = 2 * np.arcsin(mod_gap_u * 0.5 / edge_sdd)
    dalpha = dalpha_mod + dalpha_modgap

    offset = np.arctan(cfg.scanner.colOffset * colsize / sdd)
    alphas = ((np.arange(1, n_modules_x + 1) - (n_modules_x + 1) / 2.0) * dalpha + offset)
    if cfg.scanner.modDeltaAlpha != 0:
        alphas += cfg.scanner.modDeltaAlpha * dalpha

    # Module coords
    sinalphas = np.sin(alphas)
    cosalphas = np.cos(alphas)

    modcoords_x = (sdd - mod_offset_w) * sinalphas
    modcoords_y = -(sdd - mod_offset_w) * cosalphas + sid

    modwidth_z = n_rows * rowsize
    modcoords_z = ((np.arange(1, n_modules_z + 1) - (n_modules_z + 1) / 2.0) * (modwidth_z + mod_gap_v))
    voffset = np.array(cfg.scanner.modRowOffset) * rowsize
    modcoords_z = np.tile(modcoords_z[:, np.newaxis], (1, int(n_modules_x))) + np.tile(voffset, (int(n_modules_z), 1))

    modcoords_x = np.tile(modcoords_x, int(n_modules_z))
    modcoords_y = np.tile(modcoords_y, int(n_modules_z))
    modcoords_z = modcoords_z.flatten()

    modcoords = np.array([modcoords_x, modcoords_y, modcoords_z])

    # Module uvecs, vvecs
    uvecs = np.array([cosalphas, sinalphas, np.zeros(int(n_modules_x))])
    uvecs = np.tile(uvecs, (int(n_modules_z), 1))
    uvecs = uvecs.reshape(3, n_modules)

    sinbetas = np.sin(np.radians(mod_rotation_u))
    cosbetas = np.cos(np.radians(mod_rotation_u))
    sinbetas = np.tile(sinbetas, int(n_modules_x * n_modules_z))
    cosbetas = np.tile(cosbetas, int(n_modules_x * n_modules_z))
    sinalphas = np.tile(sinalphas, int(n_modules_z))
    cosalphas = np.tile(cosalphas, int(n_modules_z))
    vvecs_x = -sinbetas * sinalphas
    vvecs_y = sinbetas * cosalphas
    vvecs_z = cosbetas
    vvecs = np.array([vvecs_x, vvecs_y, vvecs_z])

    startindex = 0
    startindices = []
    n_cells = int(n_cols * n_rows)
    for i in range(n_modules):
        startindices.append(startindex)
        startindex += n_cells

    # DETECTOR DEFINITION
    if not cfg.det:
        cfg.det = CFG()

    cfg.det.n_modules = n_modules
    cfg.det.modcoords = np.single(modcoords)
    cfg.det.uvecs = np.single(uvecs)
    cfg.det.vvecs = np.single(vvecs)
    cfg.det.total_n_cells = int(n_cols * n_rows * n_modules)
    cfg.det.startindices = np.int32(startindices)
    cfg.det.n_moddefs = 1
    cfg.det.modtypes = np.ones(n_modules, dtype=np.int32)

    # MODULE DEFINITION
    # Cell and oversample parameters
    det = detector_pack_module(cfg, cfg.det)

    print("\n... done computing the detector.\n")

    return cfg.det