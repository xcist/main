# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

from gecatsim.pyfiles.GetMu import GetMu
from gecatsim.pyfiles.randpf import randpf
from gecatsim.pyfiles.CommonTools import *

def Detection_EI(cfg, viewId, subViewId):

    Evec = cfg.sim.Evec
        
    # detection efficiency
    if viewId == cfg.sim.startViewId and subViewId == 0:
        # detector prefilter
        Wvec = feval(cfg.physics.prefilterCallback, cfg)
        
        # detector absorption
        detectorMu = GetMu(cfg.scanner.detectorMaterial, Evec)
        detEff = 1-np.exp(-0.1*cfg.scanner.detectorDepth/cfg.det.cosBetas*detectorMu)
        
        cfg.sim.Wvec = Wvec*detEff
        
    # Apply energy-dependent detection efficiency
    thisSubView = cfg.thisSubView*cfg.sim.Wvec
    
    # scatter cross-talk
    if cfg.physics.crosstalkCallback:
        thisSubView = feval(cfg.physics.crosstalkCallback, thisSubView, cfg)
        
    # quantum noise
    if cfg.sim.enableQuantumNoise:
        thisSubView = randpf(thisSubView)
        
    # merge energies
    thisSubView = thisSubView.dot(Evec)
        
    # lag
    if cfg.physics.lagCallback:
        thisSubView = feval(cfg.physics.lagCallback, thisSubView, viewId, subViewId, cfg)

    # accumulate subviews
    if subViewId == 0:
        cfg.thisView = thisSubView
    else:
        cfg.thisView += thisSubView
        
    # for final subview
    if subViewId == cfg.sim.subViewCount-1:
        # optical cross-talk
        if cfg.physics.opticalCrosstalkCallback:
            cfg.thisView = feval(cfg.physics.opticalCrosstalkCallback, cfg.thisView, cfg)
        
        # DAS
        cfg.thisView = feval(cfg.physics.DASCallback, cfg.thisView, viewId, cfg)
    
    return cfg


if __name__ == "__main__":

    cfg = source_cfg("./cfg/default.cfg")
    
    cfg.sim.startViewId = 0
    cfg.sim.stopViewId = 2
    cfg.sim.enableQuantumNoise = 1
    cfg.sim.subViewCount = 1
    cfg.det.totalNumCells = 5
    cfg.det.cosBetas = np.ones([cfg.det.totalNumCells, 1], dtype=np.single)
    
    cfg = feval(cfg.protocol.spectrumCallback, cfg)
    
    cfg.thisSubView = np.float32(np.random.random([cfg.det.totalNumCells, cfg.sim.Evec.size])*100)
    cfg.thisSubView[1:4, 3:5] = 0
    viewId = 0
    subViewId = 0
    
    check_value(cfg.sim.Evec)
    check_value(cfg.thisSubView)
    
    cfg = Detection_EI(cfg, viewId, subViewId)
    check_value(cfg.sim.Wvec)
    check_value(cfg.thisView)
    
    
