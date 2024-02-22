import sys
from copy import copy, deepcopy
from gecatsim.pyfiles.CommonTools import *

def CopyCfgPhantom(cfgfrom, cfgto=None):
    if cfgto is None:
        cfgto = CFG()
    if hasattr(cfgfrom.phantom, "centerOffset"):
        cfgto.phantom.centerOffset = deepcopy(cfgfrom.phantom.centerOffset)
    if hasattr(cfgfrom.phantom, "callback"):
        cfgto.phantom.callback = deepcopy(cfgfrom.phantom.callback)
    if hasattr(cfgfrom.phantom, "filename"):
        cfgto.phantom.filename = deepcopy(cfgfrom.phantom.filename)
    if hasattr(cfgfrom.phantom, "projectorCallback"):
        cfgto.phantom.projectorCallback = deepcopy(cfgfrom.phantom.projectorCallback)
    if hasattr(cfgfrom.phantom, "projectorNumThreads"):
        cfgto.phantom.projectorNumThreads = deepcopy(cfgfrom.phantom.projectorNumThreads)
    if hasattr(cfgfrom.phantom, "scale"):
        cfgto.phantom.scale = deepcopy(cfgfrom.phantom.scale)
    if hasattr(cfgfrom.phantom, "numberOfMaterials"):
        cfgto.phantom.numberOfMaterials = deepcopy(cfgfrom.phantom.numberOfMaterials)

    return cfgto

# split cfg.phantom from a list into 
def SplitCfgPhantom(cfg):
    cfg_list = []
    for i in range(len(cfg.phantom.filename)):
        thiscfg = CFG()
        if hasattr(cfg.phantom, "centerOffset"):
            thiscfg.phantom.centerOffset = deepcopy(cfg.phantom.centerOffset[i])
        if hasattr(cfg.phantom, "callback"):
            thiscfg.phantom.callback = deepcopy(cfg.phantom.callback[i])
        if hasattr(cfg.phantom, "filename"):
            thiscfg.phantom.filename = deepcopy(cfg.phantom.filename[i])
        if hasattr(cfg.phantom, "projectorCallback"):
            thiscfg.phantom.projectorCallback = deepcopy(cfg.phantom.projectorCallback[i])
        if hasattr(cfg.phantom, "projectorNumThreads"):
            thiscfg.phantom.projectorNumThreads = deepcopy(cfg.phantom.projectorNumThreads[i])
        if hasattr(cfg.phantom, "numberOfMaterials"):
            thiscfg.phantom.numberOfMaterials = deepcopy(cfg.phantom.numberOfMaterials[i])
        if hasattr(cfg.phantom, "scale"):
            thiscfg.phantom.scale = deepcopy(cfg.phantom.scale[i])
        if hasattr(cfg.phantom, "projectorNumThreads"):
            thiscfg.phantom.projectorNumThreads = deepcopy(cfg.phantom.projectorNumThreads[i])

        cfg_list.append(thiscfg)

    return cfg_list

def PhantomWrapper(cfg):
    if isinstance(cfg.phantom.callback, list):
        # if duplicate phantom callback, raise an error
        if len(cfg.phantom.callback) != len(set(cfg.phantom.callback)):
            print("Error! XCIST does not support repeated phantom types.\n")
            sys.exit(1)

        cfg_list = SplitCfgPhantom(cfg)
        origCfgPhantom = CopyCfgPhantom(cfg)
        cfg_nom_list = []
        for thiscfg in cfg_list:
            cfg = CopyCfgPhantom(thiscfg, cfg)
            cfg = feval(cfg.phantom.callback, cfg)
            if hasattr(cfg.phantom, 'numberOfMaterials'):
                cfg_nom_list.append(cfg.phantom.numberOfMaterials)
                delattr(cfg.phantom, 'numberOfMaterials')
            else:
                # ncat does not have such property
                cfg_nom_list.append(None)

        # restore original CFG
        cfg = CopyCfgPhantom(origCfgPhantom, cfg)
        cfg.phantom.numberOfMaterials = cfg_nom_list
    else:
        cfg = feval(cfg.phantom.callback, cfg)

    return cfg

def ProjectorWrapper(cfg, viewId, subViewId):
    if isinstance(cfg.phantom.projectorCallback, list):
        cfg_list = SplitCfgPhantom(cfg)
        origCfgPhantom = CopyCfgPhantom(cfg)
        for thiscfg in cfg_list:
            cfg = CopyCfgPhantom(thiscfg, cfg)
            cfg = feval(cfg.phantom.projectorCallback, cfg, viewId, subViewId)
        # restore original CFG
        cfg = CopyCfgPhantom(origCfgPhantom, cfg)
    else:
        cfg = feval(cfg.phantom.projectorCallback, cfg, viewId, subViewId)

    return cfg
