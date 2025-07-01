# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

from gecatsim.pyfiles.CommonTools import *
from gecatsim.pyfiles.PhantomProjectorWrapper import PhantomWrapper, ProjectorWrapper

def Phantom_callback(cfg):
    if cfg.sim.isPhantomScan:
        cfg = PhantomWrapper(cfg)
    return cfg