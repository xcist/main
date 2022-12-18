# Copyright 2022, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE
import numpy as np

# only the previous view matters
def Detection_Lag(thisView, viewId, cfg):
    tau1, tau2 = cfg.physics.lag_taus
    alpha1, alpha2 = cfg.physics.lag_alphas
    dt = 1000.*cfg.protocol.rotationTime/cfg.protocol.viewsPerRotation # in ms

    if viewId == cfg.protocol.startViewId:
        cfg.prev_view = np.copy(thisView)
        return thisView

    # use two terms
    scale1 = alpha1*np.exp(-dt/tau1)/tau1
    scale2 = alpha2*np.exp(-dt/tau2)/tau2
    #print(scale1+scale2);breakpoint()
    outview = (thisView+scale1*cfg.prev_view+scale2*cfg.prev_view)/(1.+scale1+scale2)
    cfg.prev_view = np.copy(outview)

    return outview
