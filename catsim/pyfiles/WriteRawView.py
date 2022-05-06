# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

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
    dims = [cfg.scanner.detectorColCount, cfg.scanner.detectorRowCount]
    thisView = cfg.thisView.reshape(dims).T.ravel()
    
    # access mode
    if viewId == cfg.sim.startViewId:
        accessMode = 'wb'
    else:
        accessMode = 'ab'
    
    # write or append raw data
    with open(fname, accessMode) as f:
        f.write(thisView)
    
    return cfg
