# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

"""
This function reads a view_weighting.m file and assigns the view weights
based on the input collimation/pitch value.
Default has been provided in case the view_Weighting file does not exist
in the path
"""

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

    # Ensure all keys are present with default values if not found
    vw.setdefault('zsf', 1.0)
    vw.setdefault('kw', 0.0)
    vw.setdefault('beta_0', 0.825)
    vw.setdefault('beta_t', 0.175)
    vw.setdefault('vct_k', 0.0)
    vw.setdefault('vct_q', 0.0)
    vw.setdefault('vct_r', 0.0)

    return vw