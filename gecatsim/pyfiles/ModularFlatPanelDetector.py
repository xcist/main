# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import numpy as np

def modular_flat_panel_detector(cfg):
    # Shortcuts
    n_rows = cfg['rows_per_mod']
    colsize = cfg['col_size']
    rowsize = cfg['row_size']
    sid = cfg['sid']
    sdd = cfg['sdd']
    n_cols = cfg['cols_per_mod']
    n_modules = int(np.ceil(cfg['col_count'] / n_cols))
    cols_to_chop = n_modules * n_cols - cfg['col_count']

    # Cell coordinates
    cols = (np.arange(1, n_cols + 1) - (n_cols + 1) / 2.0) * colsize
    rows = (np.arange(1, n_rows + 1) - (n_rows + 1) / 2.0) * rowsize
    cellcoords = np.array([np.tile(cols, n_rows), np.repeat(rows, n_cols)])
    n_cells = n_cols * n_rows
    if cols_to_chop:
        cols = (np.arange(1, n_cols - cols_to_chop + 1) - (n_cols - cols_to_chop + 1) / 2.0) * colsize
        cellcoords = np.concatenate((cellcoords, np.array([np.tile(cols, n_rows), np.repeat(rows, n_cols - cols_to_chop)])), axis=1)
        n_cells = [n_cells, n_rows * (n_cols - cols_to_chop)]

    # Sample coordinates
    n_u = cfg['col_oversample']
    du = colsize * cfg['col_fillfraction'] / n_u
    us = (np.arange(1, n_u + 1) - (n_u + 1) / 2.0) * du
    uweights = np.ones(n_u) / n_u
    xtalk = cfg['col_crosstalk']
    if xtalk != 0:
        uweights = np.concatenate(([xtalk], uweights * (1 - 2 * xtalk), [xtalk]))
        us = np.concatenate(([us[-1] - colsize], us, [us[0] + colsize]))
        n_u += 2

    n_v = cfg['row_oversample']
    dv = rowsize * cfg['row_fillfraction'] / n_v
    vs = (np.arange(1, n_v + 1) - (n_v + 1) / 2.0) * dv
    vweights = np.ones(n_v) / n_v
    xtalk = cfg['row_crosstalk']
    if xtalk != 0:
        vweights = np.concatenate(([xtalk], vweights * (1 - 2 * xtalk), [xtalk]))
        vs = np.concatenate(([vs[-1] - rowsize], vs, [vs[0] + rowsize]))
        n_v += 2

    n_samples = n_u * n_v
    samplecoords = np.zeros((2, n_samples))
    samplecoords[0, :] = np.tile(us, n_v)
    samplecoords[1, :] = np.repeat(vs, n_u)
    weights = np.outer(vweights, uweights)

    # Module offsets
    modwidth = n_cols * colsize
    uoffset = cfg['col_offset'] * colsize
    voffset = cfg['row_offset'] * rowsize

    # Module coordinates, uvecs, vvecs
    modcoords = np.array([np.arange(-(n_modules - 1) * modwidth / 2, (n_modules - 1) * modwidth / 2 + 1, modwidth) + uoffset,
                          sid - sdd * np.ones(n_modules),
                          np.full(n_modules, voffset)])
    uvecs = np.array([np.ones(n_modules), np.zeros(n_modules), np.zeros(n_modules)])
    vvecs = np.array([np.zeros(n_modules), np.zeros(n_modules), np.ones(n_modules)])
    startindex = 0
    startindices = []
    for i in range(n_modules):
        startindices.append(startindex)
        startindex += n_cells[0] if isinstance(n_cells, list) else n_cells

    # Detector definition
    det = {
        'n_modules': n_modules,
        'modcoords': modcoords,
        'uvecs': uvecs,
        'vvecs': vvecs,
        'total_n_cells': cfg['col_count'] * cfg['rows_per_mod'],
        'startindices': startindices,
        'n_moddefs': 1 if cols_to_chop == 0 else 2,
        'modtypes': np.ones(n_modules) if cols_to_chop == 0 else np.concatenate((np.ones(n_modules - 1), [2])),
        'n_cells': n_cells,
        'cellcoords': cellcoords,
        'n_samples': np.full(1 if cols_to_chop == 0 else 2, n_samples),
        'samplecoords': np.tile(samplecoords, (1, 1, 1 if cols_to_chop == 0 else 2)),
        'weights': np.tile(weights, (1, 1, 1 if cols_to_chop == 0 else 2)),
        'activearea': np.full(1 if cols_to_chop == 0 else 2, colsize * cfg['col_fillfraction'] * rowsize * cfg['row_fillfraction']),
        'width': np.full(1 if cols_to_chop == 0 else 2, (n_cols + 1) * colsize),
        'height': np.full(1 if cols_to_chop == 0 else 2, (n_rows + 1) * rowsize)
    }

    return det
