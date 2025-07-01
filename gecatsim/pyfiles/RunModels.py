# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

from gecatsim.pyfiles.CommonTools import *
from gecatsim.pyfiles.C_Projector_SetData import C_Projector_SetData

def RunModels(cfg, viewId, subViewId, modelQueue):
    # 'detector','source','gantry','ray_angles','spectrum','filters','flux','phantom'

    if not hasattr(cfg.physics,'recalcGantry'):
        cfg.physics.recalcGantry = 2

    callback_map = {
        'detector': cfg.scanner.detectorCallback,
        'source': cfg.scanner.focalspotCallback,
        'gantry': cfg.protocol.scanTrajectory,
        'ray_angles': cfg.physics.rayAngleCallback,
        'spectrum': cfg.protocol.spectrumCallback,
        'filters': cfg.protocol.filterCallback,
        'flux': cfg.physics.fluxCallback,
        'phantom': 'Phantom_callback'
    }

    recalc_map = {
        'detector': cfg.physics.recalcDet,
        'source': cfg.physics.recalcSrc,
        'gantry': cfg.physics.recalcGantry,
        'ray_angles': cfg.physics.recalcRayAngle,
        'spectrum': cfg.physics.recalcSpec,
        'filters': cfg.physics.recalcFilt,
        'flux': cfg.physics.recalcFlux,
        'phantom': cfg.physics.recalcPht
    }

    cfg.viewId = viewId
    cfg.subViewId = subViewId
    
    for the_model in modelQueue:
        if (viewId==cfg.sim.startViewId and subViewId==0) \
            or (recalc_map[the_model]==1 and subViewId==0) \
            or recalc_map[the_model]==2:
            cfg = feval(callback_map[the_model], cfg)

    # Pass det and src to C projectors
    C_Projector_SetData(cfg, viewId, subViewId)

    return cfg
