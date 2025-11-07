# Copyright 2025, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import gecatsim as xc

def catdkvp(configfilename, cfg, adjust, preadjust, silent):
    """
    Aim
    Simulate a rotate-rotate dual kVp CT scan and write out projection data.
    """
    if hasattr(cfg, 'rotate_rotate') and cfg.rotate_rotate == 1:
        print('\n Simulating low kVp scan...\n')
        basename = cfg.results_basename + '_lowkvp'
        xc.CatSim(configfilename, cfg, f'cfg["rr_kvp"]="low_kvp"; cfg["results_basename"]="{basename}";', preadjust, silent)
        print('\n... done simulating low kVp scan.\n')

        print('\n Simulating high kVp scan...\n')
        basename = cfg.results_basename + '_highkvp'
        xc.CatSim(configfilename, cfg, f'cfg["rr_kvp"]="high_kvp"; cfg["results_basename"]="{basename}";', preadjust, silent)
        print('\n... done simulating high kVp scan.\n')