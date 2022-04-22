# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

import numpy as np
import numpy.matlib as nm
from catsim.CommonTools import *

def prep_view(cfg):
    
    ###--------- read scan data
    if hasattr(cfg.det, "totalNumCells"):
        totalNumCells = cfg.det.totalNumCells
    else:
        totalNumCells = cfg.scanner.detectorColCount*cfg.scanner.detectorRowCount
    airscan = rawread(cfg.resultsName+'.air', [cfg.protocol.airViewCount, totalNumCells], 'float')
    offsetScan = rawread(cfg.resultsName+'.offset', [cfg.protocol.offsetViewCount, totalNumCells], 'float')
    phantomScan = rawread(cfg.resultsName+'.scan', [cfg.protocol.viewCount, totalNumCells], 'float')

    ###--------- pre-log
    
    ###--------- log
    if cfg.protocol.airViewCount==1:
        airscan = nm.repmat(airscan, cfg.protocol.viewCount, 1)
    if cfg.protocol.offsetViewCount==1:
        offsetScan = nm.repmat(offsetScan, cfg.protocol.viewCount, 1)
    prep = (phantomScan-offsetScan)/(airscan-offsetScan)
    smallValue = 6.1442124e-06  # exp(-12) 
    prep[prep<smallValue] = smallValue # limits mu values to 12
    prep = -np.log(prep)
    
    ###--------- post-log

    # now perform BHC
    if hasattr(cfg.physics, "callback_post_log") and cfg.physics.callback_post_log is not None:
        print("Applying Beam Hardening Correction (ACCURATE BHC)...\n")
        prep = feval(cfg.physics.callback_post_log, cfg, prep)
        print("... done applying water BHC.\n")

    # a simple low signal correction, further limiting mu values if desired
    if cfg.protocol.maxPrep>0:
        prep[prep>cfg.protocol.maxPrep] = cfg.protocol.maxPrep
        
    ###--------- save prep
    fname = cfg.resultsName + '.prep'
    rawwrite(fname, prep)

    return cfg
