
import numpy as np
import os
from Phantom_Analytic import parse_analytical_ppm
from gecatsim.pyfiles.CommonTools import feval


def phantom_polygonal(cfg):
    print('Starting to pass POLYGONAL phantom to C...')

    if os.path.exists(cfg.phantom_filename):
        exec(open(cfg.phantom_filename).read())

    object = parse_analytical_ppm(cfg.phantom_filename)

    numPoly = len(object['type'])
    feval('C_Phantom_Polygonal_Clear', numPoly)
    
    for index in range(numPoly):
        # rotate around y-axis
        if 'phantom_rotation_y' in cfg and cfg['phantom_rotation_y'] != 0:
            ang = cfg['phantom_rotation_y']
            sina = np.sin(-ang / 180 * np.pi)
            cosa = np.cos(-ang / 180 * np.pi)
            transmat = np.array([[cosa, 0, sina], [0, 1, 0], [-sina, 0, cosa]])
            object['vertices'][index] = np.dot(object['vertices'][index], transmat)

        # Scale and position
        object['vertices'][index] = object['vertices'][index] * cfg['phantom_scale'] + np.tile(cfg['phantom_position_offset'], (object['vertices'][index].shape[0], 1))

        # Pass phantom info to C
        num_triangles = object['tri_inds'][index].shape[0]
        inds = object['tri_inds'][index].T
        triangles = object['vertices'][index][inds.ravel(), :].T.reshape((9, num_triangles))
        feval('C_Phantom_Polygonal_SetPolygon', triangles, triangles.shape[1], object['density'][index], object['material'][index] - 1)

    print('... done with phantom.')

