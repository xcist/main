# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import numpy as np

def vectornorm(xyz):
    if xyz.shape[0] != 3:
        print('ERROR: argument of vectornorm has to be of size 3 x n')
        return None

    norms = np.sqrt(np.sum(xyz * xyz, axis=0))
    return norms


def readViewWeighting(collimation, pitch):
    vw = {}
    try:
        with open('view_weighting.m', 'r') as file:
            for line in file:
                if line.strip():
                    parts = line.split()
                    vw['collimation'] = float(parts[0])
                    vw['pitch'] = float(parts[2])
                    if vw['collimation'] == collimation and vw['pitch'] == pitch:
                        vw['zsf'] = float(parts[3])
                        vw['kw'] = float(parts[4])
                        vw['beta_0'] = float(parts[5])
                        vw['beta_t'] = float(parts[6])
                        vw['vct_k'] = float(parts[7])
                        vw['vct_q'] = float(parts[8])
                        vw['vct_r'] = float(parts[9])
                        break
    except FileNotFoundError:
        # Default values if the file does not exist
        vw['beta_0'] = 0.825
        vw['beta_t'] = 0.175
        vw['vct_k'] = 0.0
        vw['vct_q'] = 0.0
        vw['vct_r'] = 0.0
        vw['zsf'] = 1.0
        vw['kw'] = 0.0

    return vw