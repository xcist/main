# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import numpy.matlib as nm
from gecatsim.pyfiles.GetMu import GetMu
from gecatsim.pyfiles.CommonTools import *

def Detection_prefilter( cfg ):
    # Wvec dim: [pixel, Ebin] ([col, row, Ebin])
    
    Evec = cfg.sim.Evec
    Wvec = np.ones(Evec.shape, dtype=np.single)
    
    if hasattr(cfg.scanner, "detectorPrefilter"):
        for ii in range(round(len(cfg.scanner.detectorPrefilter)/2)):
            material = cfg.scanner.detectorPrefilter[2*ii]
            depth = cfg.scanner.detectorPrefilter[2*ii+1]
            mu = GetMu(material, Evec)
            Wvec = Wvec*np.exp(-mu*depth*0.1)
    Wvec = nm.repmat(Wvec, cfg.det.totalNumCells, 1)
    Wvec = Wvec.astype(np.single)
    
    return Wvec
    
# if __name__ == "__main__":
#     cfg = CFG()
#     cfg.detector = CFG()
#     cfg.scanner.detectorPrefilter = ['al', 0.1, 'water', 2]
#
#     cfg = source_cfg("./cfg/default.cfg", cfg)
#
#     cfg.det.totalNumCells = 3
#
#     cfg = feval(cfg.protocol.spectrumCallback, cfg)
#     Wvec = Detection_prefilter( cfg )
#     check_value(Wvec)
