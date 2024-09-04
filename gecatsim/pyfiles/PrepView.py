# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

# This is raw-data prep for EICT or gray-scale CT only. Use other code to prep spectral CT / PCCT.

import numpy.matlib as nm
from gecatsim.pyfiles.CommonTools import *
from gecatsim.pyfiles.LowSignalCorr import LowSignalCorr

def prep_view(cfg):
    
    ###--------- read scan data
    if hasattr(cfg.det, "totalNumCells"):
        totalNumCells = cfg.det.totalNumCells
    else:
        totalNumCells = cfg.scanner.detectorColCount*cfg.scanner.detectorRowCount
    airscan = rawread(cfg.resultsName+'.air', [cfg.protocol.airViewCount, totalNumCells], 'float')
    offsetScan = rawread(cfg.resultsName+'.offset', [cfg.protocol.offsetViewCount, totalNumCells], 'float')
    phantomScan = rawread(cfg.resultsName+'.scan', [cfg.protocol.viewCount, totalNumCells], 'float')

    if cfg.protocol.airViewCount==1:
        airscan = nm.repmat(airscan, cfg.protocol.viewCount, 1)
    if cfg.protocol.offsetViewCount==1:
        offsetScan = nm.repmat(offsetScan, cfg.protocol.viewCount, 1)
        
    ###--------- pre-log
    if hasattr(cfg.physics, "callback_pre_log") and cfg.physics.callback_pre_log:
        airscan, offsetScan, phantomScan = feval(cfg.physics.callback_pre_log, cfg, airscan, offsetScan, phantomScan)
    
    
    ###--------- log
    prep = (phantomScan-offsetScan)/(airscan-offsetScan)
    
    ### simple low-signal correction and -log
    #smallValue = 1e-12
    #prep[prep<smallValue] = smallValue
    #prep = -np.log(prep)
    
    ### signal-domain low-signal correction and -log
    prep = LowSignalCorr(cfg, prep)
    prep[prep<0] = 0
    
    
    ###--------- post-log
    # BHC and so furth
    if hasattr(cfg.physics, "callback_post_log") and cfg.physics.callback_post_log:
        prep = feval(cfg.physics.callback_post_log, cfg, prep)
        
    # a simple p-domain low-signal correction, further limiting mu values if desired
    if cfg.protocol.maxPrep>0:
        prep[prep>cfg.protocol.maxPrep] = cfg.protocol.maxPrep
        
    ###--------- save prep
    fname = cfg.resultsName + '.prep'
    rawwrite(fname, prep)

    return cfg
