# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

"""
Aim
    Returns a modular flat panel detector. Only validated for ceil(cfg.col_count/cfg.cols_per_mod) = cfg.col_count/cfg.cols_per_mod. Modules make it possible to use parallelized projection code.

Inputs
    cfg.cols_per_mod : number of detector columns per module
                       if this does not evenly divide into cfg.col_count, the last module will have fewer columns
    cfg.rows_per_mod : number of detector rows per module
    cfg.total_n_cells : total number of cells in the detector
    cfg.col_count : total number of detector columns
    cfg.col_size : physical detector column size (pitch, in mm)
    cfg.row_size : physical detector row size (pitch, in mm)
    cfg.sid : source-to-iso distance (in mm) for a source-focused detector
    cfg.sdd : source-to-detector distance (in mm) for a source-focused detector
    cfg.col_oversample : number of detector column subsamples
                         2 more samples will be added for x-crosstalk
    cfg.row_oversample : number of detector row subsamples
                         2 more samples will be added for y-crosstalk
    cfg.col_crosstalk : fractional crosstalk to each neighboring column
                        eg: 0.02 results in 96% weight for current column
                        and 2% for each neighboring column
    cfg.row_crosstalk : fractional crosstalk to each neighboring row
                        eg: 0.02 results in 96% weight for current row
                        and 2% for each neighboring row
    cfg.col_fillfraction : column fill fraction
    cfg.row_fillfraction : row fill fraction
    cfg.col_offset : in-plane detector offset (in alpha)
    cfg.row_offset : longitudinal detector offset (in z)

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
"""
import numpy as np
from gecatsim.pyfiles.CommonTools import *

def modular_flat_panel_detector(cfg):
    print("\nComputing a MODULAR FLAT PANEL detector...\n")

    # SHORTCUTS
    sid = cfg.scanner.sid
    sdd = cfg.scanner.sdd
    n_rows = cfg.scanner.detectorRowsPerMod
    colsize = cfg.scanner.detectorColSize
    rowsize = cfg.scanner.detectorRowSize
    n_cols = cfg.scanner.detectorColsPerMod
    n_modules = int(np.ceil(cfg.scanner.detectorColCount / n_cols))
    cols_to_chop = n_modules * n_cols - cfg.scanner.detectorColCount

    # CELL.COORDS
    cols = (np.arange(1, n_cols + 1) - (n_cols + 1) / 2.0) * colsize
    rows = (np.arange(1, n_rows + 1) - (n_rows + 1) / 2.0) * rowsize
    tmp_var = np.tile(cols, (n_rows, 1))
    cellcoords = np.array([tmp_var.flatten(), np.repeat(rows, n_cols)])

    n_cells = [n_cols * n_rows]

    if cols_to_chop:
        cols = (np.arange(1, n_cols - cols_to_chop + 1) - (n_cols - cols_to_chop + 1) / 2.0) * colsize
        tmp_var = np.tile(cols, (n_rows, 1))
        cellcoords = np.concatenate((cellcoords[:, :, np.newaxis],
                                     np.array([tmp_var.flatten(), np.repeat(rows, n_cols - cols_to_chop)])[:, :, np.newaxis]), axis=2)
        n_cells.append(n_rows * (n_cols - cols_to_chop))

    # SAMPLE.U COORDS
    n_u = cfg.scanner.colOversample
    du = colsize * cfg.scanner.colFillFraction / n_u
    us = (np.arange(1, n_u + 1) - (n_u + 1) / 2.0) * du
    uweights = np.ones(n_u) / n_u

    xtalk = cfg.scanner.colCrosstalk

    if xtalk != 0:
        uweights = np.concatenate(([xtalk], uweights * (1 - 2 * xtalk), [xtalk]))
        if xtalk > uweights[1]:
            print('xtalk too high for uniform sample weights')
            return None

        us = np.concatenate(([us[-1] - colsize], us, [us[0] + colsize]))
        n_u += 2

    # SAMPLE.V COORDS
    n_v = cfg.scanner.rowOversample
    dv = rowsize * cfg.scanner.rowFillFraction / n_v
    vs = (np.arange(1, n_v + 1) - (n_v + 1) / 2.0) * dv

    vweights = np.ones(n_v) / n_v

    xtalk = cfg.scanner.rowCrosstalk

    if xtalk != 0:
        vweights = np.concatenate(([xtalk], vweights * (1 - 2 * xtalk), [xtalk]))
        if xtalk > vweights[1]:
            print('xtalk too high for uniform sample weights')
            return None

        vs = np.concatenate(([vs[-1] - rowsize], vs, [vs[0] + rowsize]))
        n_v += 2

    # SAMPLE.COORDS
    n_samples = n_u * n_v
    samplecoords = np.zeros((2, n_samples))
    tmp_var = np.tile(us, (n_v, 1))
    samplecoords[0, :] = tmp_var.flatten()
    tmp_var = np.tile(vs[:, np.newaxis], (1, n_u))
    samplecoords[1, :] = tmp_var.flatten()
    weights = np.outer(vweights, uweights)

    # MODULE.OFFSETS
    modwidth = n_cols * colsize
    uoffset = cfg.scanner.colOffset * colsize
    voffset = cfg.scanner.rowOffset * rowsize

    # MODULE.COORDS, UVECS, VVECS
    modcoords = np.array(
        [np.arange(-(n_modules - 1) * modwidth / 2, (n_modules - 1) * modwidth / 2 + modwidth, modwidth) + uoffset,
         sid - sdd * np.ones(n_modules),
         np.tile(voffset, n_modules)])
    uvecs = np.array([np.ones(n_modules), np.zeros(n_modules), np.zeros(n_modules)])
    vvecs = np.array([np.zeros(n_modules), np.zeros(n_modules), np.ones(n_modules)])
    startindex = 0
    startindices = []
    for i in range(n_modules):
        startindices.append(startindex)
        startindex += n_cells[0]

    # Detector definition
    if not cfg.det:
        cfg.det = CFG()

    cfg.det.n_modules = n_modules
    cfg.det.modcoords = modcoords
    cfg.det.uvecs = uvecs
    cfg.det.vvecs = vvecs
    cfg.det.total_n_cells = cfg.scanner.detectorColCount * cfg.scanner.detectorRowsPerMod
    cfg.det.startindices = startindices
    cfg.det.n_moddefs = 1 if cols_to_chop == 0 else 2
    cfg.det.modtypes = np.ones(n_modules) if cols_to_chop == 0 else np.concatenate((np.ones(n_modules - 1), [2]))
    cfg.det.n_cells = n_cells
    cfg.det.cellcoords = cellcoords
    cfg.det.n_samples = np.full(1 if cols_to_chop == 0 else 2, n_samples)
    cfg.det.samplecoords = np.tile(samplecoords[:, :, np.newaxis], (1, 1, 1 if cols_to_chop == 0 else 2))
    cfg.det.weights = np.tile(weights[:, :, np.newaxis], (1, 1, 1 if cols_to_chop == 0 else 2))
    cfg.det.activearea = np.full(1 if cols_to_chop == 0 else 2, colsize * cfg.scanner.colFillFraction * rowsize * cfg.scanner.rowFillFraction)
    cfg.det.width = np.full(1 if cols_to_chop == 0 else 2, (n_cols + 1) * colsize)
    cfg.det.height = np.full(1 if cols_to_chop == 0 else 2, (n_rows + 1) * rowsize)

    return cfg.det