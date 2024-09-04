# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license
import numpy as np

# only the previous view matters
def Detection_Lag(thisView, viewId, subViewId, cfg):
    tau1, tau2 = cfg.physics.lag_taus
    alpha1, alpha2 = cfg.physics.lag_alphas
    dt = 1000.*cfg.protocol.rotationTime/cfg.protocol.viewsPerRotation # in ms
    dt /= cfg.physics.viewSampleCount

    unaccounted = 1 - alpha1 - alpha2
    invintegral = alpha1*(1.0 - np.exp(-dt/2/tau1)) + alpha2*(1.0 - np.exp(-dt/2/tau2)) + unaccounted

    if viewId == cfg.protocol.startViewId and subViewId == 0:
        if cfg.sim.isAirScan:
            cfg.memview1 = thisView*np.exp(-dt/2/tau1)/(1-np.exp(-dt/tau1))
            cfg.memview2 = thisView*np.exp(-dt/2/tau2)/(1-np.exp(-dt/tau2))
        #else cfg.sim.isPhantomScan:
        else:
            cfg.memview1 = 0
            cfg.memview2 = 0
        #return thisView*invintegral

    outview = invintegral*thisView + (alpha1*(1.0 - np.exp(-dt/tau1)))*cfg.memview1 + (alpha2*(1.0 - np.exp(-dt/tau2)))*cfg.memview2
    cfg.memview1 = (cfg.memview1 * np.exp(-dt/tau1) + thisView*np.exp(-dt/2/tau1))
    cfg.memview2 = (cfg.memview2 * np.exp(-dt/tau2) + thisView*np.exp(-dt/2/tau2))

    return outview
