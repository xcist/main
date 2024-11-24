# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import numpy as np
from gecatsim.pyfiles.USampling import USampling

def detector_pack_module(cfg, det):
    global Verbose

    print("\nComputing cell sample coordinates for PACK MODULES...\n")

    colsize = cfg['col_size']
    rowsize = cfg['row_size']

    n_cols = cfg['cols_per_mod']
    n_rows = cfg['rows_per_mod']

    n_packs_u = cfg['col_packs_per_mod']
    n_packs_v = cfg['row_packs_per_mod']

    pack_gap_u = cfg['col_pack_gap']
    pack_gap_v = cfg['row_pack_gap']

    # CELL.COORDS
    xcells = n_cols / n_packs_u
    zcells = n_rows / n_packs_v
    cols = np.zeros(n_cols)
    rows = np.zeros(n_rows)

    for n in range(n_cols):
        packindex = np.floor((n) / xcells)
        cols[n] = (n - (n_cols + 1) / 2) * colsize + (packindex - (n_packs_u - 1) / 2) * pack_gap_u

    for n in range(n_rows):
        packindex = np.floor((n) / zcells)
        rows[n] = (n - (n_rows + 1) / 2) * rowsize + (packindex - (n_packs_v - 1) / 2) * pack_gap_v

    tmp_rep_cols = np.tile(cols, (n_rows, 1))
    cellcoords = np.vstack((tmp_rep_cols.flatten(), np.repeat(rows, n_cols)))
    n_cells = n_cols * n_rows

    # SAMPLE.U COORDS
    n_u = cfg['col_oversample']
    if 'col_intensive_oversample' not in cfg or cfg['col_intensive_oversample'] < 0 or \
            'col_intensive_oversample_length' not in cfg or cfg['col_intensive_oversample_length'] < 0:
        cfg['col_intensive_oversample'] = 0
        cfg['col_intensive_oversample_length'] = 0

    intensive_n = cfg['col_intensive_oversample']
    us, u_step = USampling(cfg, n_u, intensive_n)
    uactive = np.sum(u_step)
    du = uactive / n_u

    if intensive_n == 0:
        uweights = np.ones(n_u) / n_u
    else:
        intensive_len = cfg['col_intensive_oversample_length']
        uweights = np.zeros(n_u)
        uweights[:intensive_n] = np.ones(intensive_n) * intensive_len / intensive_n / uactive
        uweights[-intensive_n:] = uweights[:intensive_n]
        uweights[intensive_n:-intensive_n] = np.ones(n_u - 2 * intensive_n) * (
                    1 - np.sum(uweights[:intensive_n]) * 2) / (n_u - 2 * intensive_n)

    xtalk = cfg['col_crosstalk']
    if xtalk != 0:
        uweights = np.concatenate(([xtalk], uweights * (1 - 2 * xtalk), [xtalk]))
        if xtalk > 0.5:
            print('column xtalk too high')
            return None
        us = np.concatenate(([-colsize], us, [colsize]))
        n_u += 2

    # SAMPLE.V COORDS
    n_v = cfg['row_oversample']
    vactive = rowsize - cfg['row_cast']
    dv = vactive / n_v
    vs = ((np.arange(1, n_v + 1) - (n_v + 1) / 2) * dv)
    vweights = np.ones(n_v) / n_v

    xtalk = cfg['row_crosstalk']
    if xtalk != 0:
        vweights = np.concatenate(([xtalk], vweights * (1 - 2 * xtalk), [xtalk]))
        if xtalk > 0.5:
            print('row xtalk too high')
            return None
        vs = np.concatenate(([-rowsize], vs, [rowsize]))
        n_v += 2

    # SAMPLE.COORDS and WEIGHTS
    n_samples = n_u * n_v
    samplecoords = np.zeros((2, n_samples))

    tmp_rpus = np.tile(us, (n_v, 1))
    tmp_rpvs = np.tile(vs[:, np.newaxis], (1, n_u))
    samplecoords[0, :] = tmp_rpus.flatten()
    samplecoords[1, :] = tmp_rpvs.flatten()

    weights = np.outer(vweights, uweights)

    # MODULE DEFINITION
    det['n_cells'] = n_cells
    det['cellcoords'] = cellcoords
    det['n_samples'] = n_samples
    det['samplecoords'] = samplecoords
    det['weights'] = weights
    det['activearea'] = uactive * vactive
    det['width'] = (n_cols + 1) * colsize
    det['height'] = (n_rows + 1) * rowsize

    det['cols'] = cols
    det['rows'] = rows
    det['sample_du'] = du
    det['sample_dv'] = dv

    print("\n... done computing cell sample coordinates.\n")

    return det