# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import numpy as np
import numpy.matlib as nm
from gecatsim.pyfiles.CommonTools import *

def Detection_Flux(cfg):
    '''
    Compute the photon flux at the detector per cell per subview
    Flux dim: [pixel, Ebin] ([col, row, Ebin])
    Note: the unit of spectrum: photons per mAs per mm^2 at 1-m distance
          and spec.Ivec is already scaled to mA and view time, i.e. mAs
    Mingye Wu, GE Research
    
    '''
    ###------- offset scan, flux = 0
    #if hasattr(cfg, 'sim.isOffsetScan') and cfg.sim.isOffsetScan:
    if cfg.sim.isOffsetScan:
        cfg.detFlux = np.zeros([cfg.det.totalNumCells, cfg.spec.nEbin], dtype=np.single)
        return cfg
    
    ###------- air or phantom scan
    detActiveArea = cfg.det.activeArea*np.cos(cfg.det.betas) # mm^2
    #detActiveArea = nm.repmat(detActiveArea, 1, cfg.spec.nEbin)
    
    distanceFactor = np.square(1000/cfg.det.rayDistance) # mm
    #distanceFactor = nm.repmat(distanceFactor, 1, cfg.spec.nEbin)
    
    cfg.spec.netIvec = cfg.spec.Ivec*cfg.src.filterTrans
    cfg.detFlux = np.multiply(cfg.spec.netIvec.astype(np.single, copy=False), np.single(detActiveArea*distanceFactor))
    
    return cfg

# if __name__ == "__main__":
#
#     cfg = source_cfg("./cfg/default.cfg")
#
#     cfg.sim.isOffsetScan = 0
#
#     cfg = feval(cfg.scanner.detectorCallback, cfg)
#     cfg = feval(cfg.scanner.focalspotCallback, cfg)
#     cfg = feval(cfg.physics.rayAngleCallback, cfg)
#     cfg = feval(cfg.protocol.spectrumCallback, cfg)
#     cfg = feval(cfg.protocol.filterCallback, cfg)
#
#     cfg = Detection_Flux(cfg)
#     check_value(cfg.detFlux)
