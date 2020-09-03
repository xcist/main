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
    
    # access mode
    if viewId == cfg.sim.startViewId:
        accessMode = 'wb'
    else:
        accessMode = 'ab'
    
    # write or append raw data
    with open(fname, accessMode) as f:
        f.write(cfg.thisView)
    
    return cfg
