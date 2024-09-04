# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

'''
Write view data to file (binary).
Mingye Wu, GE Research

'''
def WriteRawView(cfg, viewId):
    # filename
    scanTypeInd = [cfg.sim.isAirScan, cfg.sim.isOffsetScan, cfg.sim.isPhantomScan].index(1)
    extName = ['.air', '.offset', '.scan'][scanTypeInd]
    fname = cfg.resultsName + extName
    
    # Change the dim from row->col to col->row
    if cfg.thisView.ndim==1:
        dims = [cfg.scanner.detectorColCount, cfg.scanner.detectorRowCount]
        thisView = cfg.thisView.reshape(dims).T.ravel()
    else:
        dims = [cfg.scanner.detectorColCount, cfg.scanner.detectorRowCount, cfg.thisView.shape[1]]  # [col,row,bin]
        thisView = cfg.thisView.reshape(dims).transpose((1,0,2)).ravel()
    
    # access mode
    if viewId == cfg.sim.startViewId:
        accessMode = 'wb'
    else:
        accessMode = 'ab'
    
    # write or append raw data
    with open(fname, accessMode) as f:
        f.write(thisView)
    
    return cfg
