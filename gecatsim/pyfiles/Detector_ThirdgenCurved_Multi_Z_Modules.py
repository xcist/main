# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import numpy as np
from gecatsim.pyfiles.Detector_Pack_Module import Detector_Pack_Module

def detector_thirdgen_curved_multi_z_modules(cfg):
    global Verbose

    print("\nComputing a THIRD-GENERATION CURVED detector with MULTIPLE MODULES IN Z...\n")

    # SHORTCUTS
    sid = float(cfg['sid'])
    sdd = float(cfg['sdd'])

    colsize = float(cfg['col_size'])
    rowsize = float(cfg['row_size'])
    n_rows = float(cfg['rows_per_mod'])
    n_cols = float(cfg['cols_per_mod'])

    n_modules_x = float(cfg['mods_per_det_x'])
    n_modules_z = float(cfg['mods_per_det_z'])
    n_modules = int(n_modules_x * n_modules_z)

    mod_gap_u = float(cfg['mod_gap_u'])
    mod_gap_v = float(cfg['mod_gap_v'])

    mod_offset_w = float(cfg['mod_offset_w'])
    mod_rotation_u = float(cfg['mod_rotation_u'])

    # Module alphas in X-Y plane
    modwidth_x = n_cols * colsize + cfg['col_pack_gap'] * (cfg['col_packs_per_mod'] - 1)
    dalpha_mod = 2 * np.arctan(modwidth_x * 0.5 / cfg['sdd'])
    edge_sdd = np.sqrt((modwidth_x * 0.5) ** 2 + cfg['sdd'] ** 2)
    dalpha_modgap = 2 * np.arcsin(mod_gap_u * 0.5 / edge_sdd)
    dalpha = dalpha_mod + dalpha_modgap

    offset = np.arctan(cfg['col_offset'] * colsize / sdd)
    alphas = ((np.arange(1, n_modules_x + 1) - (n_modules_x + 1) / 2.0) * dalpha + offset)
    if 'mod_delta_alpha' in cfg and cfg['mod_delta_alpha'] != 0:
        alphas += cfg['mod_delta_alpha'] * dalpha

    # Module coords
    sinalphas = np.sin(alphas)
    cosalphas = np.cos(alphas)

    modcoords_x = (sdd - mod_offset_w) * sinalphas
    modcoords_y = -(sdd - mod_offset_w) * cosalphas + sid

    modwidth_z = n_rows * rowsize
    modcoords_z = ((np.arange(1, n_modules_z + 1) - (n_modules_z + 1) / 2.0) * (modwidth_z + mod_gap_v))
    voffset = cfg['mod_row_offset'] * rowsize
    modcoords_z = np.tile(modcoords_z[:, np.newaxis], (1, int(n_modules_x))) + np.tile(voffset, (int(n_modules_z), 1))

    modcoords = np.vstack((modcoords_x.flatten(), modcoords_y.flatten(), modcoords_z.flatten()))

    # Module uvecs, vvecs
    uvecs = np.vstack((cosalphas, sinalphas, np.zeros(int(n_modules_x))))
    uvecs = np.tile(uvecs, (int(n_modules_z), 1)).reshape(3, n_modules)

    sinbetas = np.sin(mod_rotation_u / 180 * np.pi)
    cosbetas = np.cos(mod_rotation_u / 180 * np.pi)
    vvecs_x = -sinbetas[:, np.newaxis] * sinalphas
    vvecs_y = sinbetas[:, np.newaxis] * cosalphas
    vvecs_z = np.tile(cosbetas[:, np.newaxis], (1, int(n_modules_x)))
    vvecs = np.vstack((vvecs_x.flatten(), vvecs_y.flatten(), vvecs_z.flatten()))

    startindex = 0
    startindices = []
    n_cells = int(cfg['cols_per_mod'] * cfg['rows_per_mod'])
    for i in range(n_modules):
        startindices.append(startindex)
        startindex += n_cells

    # DETECTOR DEFINITION
    det = {}

    det['n_modules'] = n_modules
    det['modcoords'] = modcoords
    det['uvecs'] = uvecs
    det['vvecs'] = vvecs
    det['total_n_cells'] = int(n_cols * n_rows * n_modules)
    det['startindices'] = startindices
    det['n_moddefs'] = 1
    det['modtypes'] = np.ones(n_modules, dtype=int)

    # MODULE DEFINITION
    # Cell and oversample parameters

    det = detector_pack_module(cfg, det)

    print("\n... done computing the detector.\n")

    return det