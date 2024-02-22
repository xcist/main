# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

import copy, time
import matplotlib.pyplot as plt
from gecatsim.pyfiles.CommonTools import *
from tqdm import tqdm
from gecatsim.pyfiles.PhantomProjectorWrapper import PhantomWrapper, ProjectorWrapper

def one_scan(cfg):
    cfg = initialize_scan(cfg)
    
    # view loop
    if cfg.sim.isPhantomScan:
        print('phantom scan view loop...')
        tmp = tqdm(range(cfg.sim.startViewId, cfg.sim.stopViewId+1))
    else:
        tmp = range(cfg.sim.startViewId, cfg.sim.stopViewId+1)
    for viewId in tmp:
    #for viewId in range(cfg.sim.startViewId, cfg.sim.stopViewId+1):
        # detector
        if viewId == cfg.sim.startViewId or cfg.physics.recalcDet:
            cfg = feval(cfg.scanner.detectorCallback, cfg)

        # source
        if viewId == cfg.sim.startViewId or cfg.physics.recalcSrc:
            cfg = feval(cfg.scanner.focalspotCallback, cfg)

        # ray angles
        if viewId == cfg.sim.startViewId or cfg.physics.recalcRayAngle:
            cfg = feval(cfg.physics.rayAngleCallback, cfg)

        # spectrum
        if viewId == cfg.sim.startViewId or cfg.physics.recalcSpec:
            cfg = feval(cfg.protocol.spectrumCallback, cfg)

        # filters (bowtie and flat)
        if viewId == cfg.sim.startViewId or cfg.physics.recalcFilt:
            cfg = feval(cfg.protocol.filterCallback, cfg)

        # flux
        if viewId == cfg.sim.startViewId or cfg.physics.recalcFlux:
            cfg = feval(cfg.physics.fluxCallback, cfg)

        # phantom and material
        if (viewId == cfg.sim.startViewId or cfg.physics.recalcPht) and cfg.sim.isPhantomScan:
            cfg = PhantomWrapper(cfg)
    
        for subViewId in range(cfg.sim.subViewCount):
            # initial subview
            cfg.thisSubView = copy.copy(cfg.detFlux)
            
            if cfg.sim.isPhantomScan:
                # apply gantry rotation and table movement
                cfg = feval(cfg.protocol.scanTrajectory, cfg, viewId)

                # projector
                cfg = ProjectorWrapper(cfg, viewId, subViewId)

                # scatter
                if cfg.physics.scatterCallback:
                    cfg = feval(cfg.physics.scatterCallback, cfg, viewId, subViewId)
            
            # detection, the input: cfg.thisSubView [col-row, energy]
            #            the output: cfg.thisView [col, row]
            cfg = feval(cfg.scanner.detectionCallback, cfg, viewId, subViewId)

            cfg = update_scan_time(cfg, subViewId)

        # save cfg.thisView to file
        cfg = feval(cfg.physics.outputCallback, cfg, viewId)
        
    #print("Scan sim time: %.1f s" % (time.time()-cfg.sim.timer))
    
    return cfg

def initialize_scan(cfg):
    # scan type
    cfg.sim.isAirScan = cfg.sim.thisScanType[0]
    cfg.sim.isOffsetScan = cfg.sim.thisScanType[1]
    cfg.sim.isPhantomScan = not(cfg.sim.isAirScan or cfg.sim.isOffsetScan)
    
    if cfg.sim.isPhantomScan:
        cfg.sim.startViewId = cfg.protocol.startViewId
        cfg.sim.stopViewId = cfg.protocol.stopViewId
        cfg.sim.subViewCount = cfg.physics.viewSampleCount
        cfg.sim.enableQuantumNoise = cfg.physics.enableQuantumNoise
        if cfg.physics.enableElectronicNoise:
            cfg.sim.eNoise = cfg.scanner.eNoise
        else:
            cfg.sim.eNoise = 0
    else:
        cfg.sim.startViewId = 0
        if cfg.sim.isAirScan:
            cfg.sim.stopViewId = cfg.protocol.airViewCount-1
        else:
            cfg.sim.stopViewId = cfg.protocol.offsetViewCount-1
        cfg.sim.subViewCount = 1
        cfg.sim.enableQuantumNoise = 0
        cfg.sim.eNoise = 0
    cfg.sim.viewCount = cfg.sim.stopViewId-cfg.sim.startViewId+1
    
    # scan time
    startTime = 0
    viewTime = cfg.protocol.rotationTime/cfg.protocol.viewsPerRotation  # view interval
    subViewTime = viewTime*cfg.protocol.dutyRatio/cfg.sim.subViewCount  # beam-on time oversampling
    cfg.time = startTime + cfg.protocol.startViewId*viewTime - 0.5*viewTime + 0.5*subViewTime
    cfg.viewTime = viewTime
    cfg.subViewTime = subViewTime
    
    # simulation time
    cfg.sim.timer = time.time()
    
    return cfg

def update_scan_time(cfg, subViewId):
    cfg.time += cfg.subViewTime
    if cfg.protocol.dutyRatio<1 and subViewId==cfg.sim.subViewCount-1:
        cfg.time = cfg.time + (1-cfg.protocol.dutyRatio) * cfg.viewTime

    return cfg

if __name__ == "__main__":

    cfg = source_cfg("./cfg/default.cfg")
    
    cfg.sim.thisScanType = [0, 0, 1]
    cfg.sim.enableQuantumNoise = 0
    cfg.sim.subViewCount = 3
    cfg.protocol.startViewId = 0
    cfg.protocol.stopViewId = cfg.protocol.viewsPerRotation-1
    cfg.protocol.startAngle = 29
    cfg.protocol.startZ = 2
    
    cfg = one_scan(cfg)
    
    thisScan = cfg.thisView.reshape(cfg.scanner.detectorColCount, cfg.scanner.detectorRowCount)

    check_value(thisScan[400:500, 7])
    plt.plot(thisScan[:, 7])
    plt.show()
