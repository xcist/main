# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

'''
Write view data to file (binary).
Mingye Wu, GE Research

'''
def WriteRawViewChunk(cfg, viewId):
    # access mode
    if viewId == cfg.sim.startViewId:
        accessMode = 'wb'
        cfg.dump_views = b''
    elif (viewId-cfg.sim.startViewId)%cfg.physics.dump_period==0 or viewId == cfg.sim.stopViewId:
        accessMode = 'ab'
    else:
        accessMode = None

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
    
    cfg.dump_views += thisView.tobytes()
    
    # write or append raw data
    if accessMode is None: return cfg

    with open(fname, accessMode) as f:
        f.write(cfg.dump_views)
    cfg.dump_views = b''
    
    return cfg
