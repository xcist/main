# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import numpy as np

def detector_skew_and_offset(cfg, det):
    # Individual Cell/Column Offsets
    if 'individual_cell_offset' in cfg and cfg['individual_cell_offset'] == 1:
        cell_x_offset = np.zeros(cfg['col_count'])
        cell_z_offset = np.zeros(cfg['col_count'])

        catsim_filename = 'catsim'  # Placeholder for the actual path
        basepath = catsim_filename.replace('\\', '/').replace('base/mfiles/catsim.m', '')
        cell_offset_directory = f"{basepath}SVCT/Cell_offset/"

        if 'cell_offset_data_file' in cfg and cfg['cell_offset_data_file']:
            cell_offset_filename = f"{cell_offset_directory}{cfg['cell_offset_data_file']}"
            try:
                with open(cell_offset_filename, 'r') as fp:
                    tline = fp.readline()
                    tmp = np.loadtxt(fp)
                    num_tline = np.fromstring(tline, sep=' ')
                    if num_tline.size == 2:
                        tmp = np.hstack((num_tline, tmp))
                    if tmp.size == 2 * cfg['col_count']:
                        tmp = tmp.reshape((2, cfg['col_count']))
                        cell_x_offset = tmp[0, :]
                        cell_z_offset = tmp[1, :]
                    else:
                        raise ValueError('Cell offset data file error')
            except FileNotFoundError:
                raise FileNotFoundError('Cannot find cell offset data file')

        cell_offset = np.array([cell_x_offset, cell_z_offset])
        cell_offset = np.tile(cell_offset[:, :, np.newaxis], (1, 1, cfg['row_count']))
        cell_offset = np.transpose(cell_offset, (0, 2, 1)).reshape((2, det['n_cells'], det['n_modules']))

        cellcoords_tmp = np.zeros_like(cell_offset)
        for mod_id in range(det['n_modules']):
            mod_type_id = det['modtypes'][mod_id]
            cellcoords_tmp[:, :, mod_id] = det['cellcoords'][:, :, mod_type_id] + cell_offset[:, :, mod_id]

        det['cellcoords'] = cellcoords_tmp
        det['n_moddefs'] = det['n_modules']
        det['modtypes'] = np.arange(1, det['n_modules'] + 1)
        det['n_cells'] = np.tile(det['n_cells'][0, 0], (1, det['n_moddefs']))
        det['n_samples'] = np.tile(det['n_samples'][0, 0], (1, det['n_moddefs']))
        det['samplecoords'] = np.tile(det['samplecoords'][:, :, 0], (1, 1, det['n_moddefs']))
        det['weights'] = np.tile(det['weights'][:, :, 0], (1, 1, det['n_moddefs']))
        det['activearea'] = np.tile(det['activearea'][0, 0], (1, det['n_moddefs']))
        det['width'] = np.tile(det['width'][0, 0], (1, det['n_moddefs']))
        det['height'] = np.tile(det['height'][0, 0], (1, det['n_moddefs']))

    # Individual Pack Skew
    if 'individual_pack_skew' in cfg and cfg['individual_pack_skew'] == 1:
        if 'pack_skew_w' not in cfg:
            raise ValueError('Pack Skew Angle Not Defined')
        elif len(cfg['pack_skew_w']) != det['n_modules'] * cfg['col_packs_per_mod']:
            raise ValueError('Number of Pack Skew Angles should be total packs per detector')

        pack_skew_w = np.radians(cfg['pack_skew_w'])
        cells_per_pack = cfg['cols_per_pack'] * cfg['rows_per_pack']
        cell_xz_skewed = np.zeros((2, cells_per_pack * cfg['col_packs_per_mod'], det['n_modules']))

        for mod_id in range(det['n_modules']):
            mod_type_id = det['modtypes'][mod_id]
            for pack_id in range(cfg['col_packs_per_mod']):
                skew = pack_skew_w[mod_id * cfg['col_packs_per_mod'] + pack_id]
                rot_mat = np.array([[np.cos(skew), -np.sin(skew)], [np.sin(skew), np.cos(skew)]])
                pack_xz = np.array([(-cfg['col_packs_per_mod'] / 2 - 0.5 + pack_id) * (
                            cfg['cols_per_pack'] * cfg['col_size'] + cfg['col_pack_gap']), 0])
                for cell_id in range(cells_per_pack):
                    ind = pack_id * cells_per_pack + cell_id
                    cell_xz = det['cellcoords'][:, ind, mod_type_id]
                    cell_xz_skewed[:, ind, mod_id] = rot_mat @ (cell_xz - pack_xz) + pack_xz

        det['cellcoords'] = cell_xz_skewed
        det['n_moddefs'] = det['n_modules']
        det['modtypes'] = np.arange(1, det['n_modules'] + 1)
        det['n_cells'] = np.tile(det['n_cells'][0, 0], (1, det['n_moddefs']))
        det['n_samples'] = np.tile(det['n_samples'][0, 0], (1, det['n_moddefs']))
        det['samplecoords'] = np.tile(det['samplecoords'][:, :, 0], (1, 1, det['n_moddefs']))
        det['weights'] = np.tile(det['weights'][:, :, 0], (1, 1, det['n_moddefs']))
        det['activearea'] = np.tile(det['activearea'][0, 0], (1, det['n_moddefs']))
        det['width'] = np.tile(det['width'][0, 0], (1, det['n_moddefs']))
        det['height'] = np.tile(det['height'][0, 0], (1, det['n_moddefs']))

    # Individual Module Skew
    if 'individual_module_skew' in cfg and cfg['individual_module_skew'] == 1:
        if 'module_skew_w' not in cfg:
            raise ValueError('Module Skew Angle Not Defined')
        elif len(cfg['module_skew_w']) != det['n_modules']:
            raise ValueError('Number of Module Skew Angles should be total modules per detector')

        mod_skew_w = np.radians(cfg['module_skew_w'])
        uv_0 = det['uvecs']
        vv_0 = det['vvecs']

        uv_1 = np.vstack((uv_0[:2, :] * np.cos(mod_skew_w), np.sin(mod_skew_w)))
        vv_1 = np.vstack((-uv_0[:2, :] * np.sin(mod_skew_w), np.cos(mod_skew_w)))

        det['uvecs'] = uv_1
        det['vvecs'] = vv_1

    return det