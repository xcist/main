import sys
from copy import copy, deepcopy
from gecatsim.pyfiles.CommonTools import *

def CopyCfgPhantom(cfgfrom, cfgto=None):
    if cfgto is None:
        cfgto = CFG()
    cfgto.phantom = deepcopy(cfgfrom.phantom)

    return cfgto

# split cfg.phantom from a list into 
def SplitCfgPhantom(cfg):
    cfg_list = []
    cfgphantom_attrs = [x for x in dir(cfg.phantom) if not x.startswith('__')]
    for i in range(len(cfg.phantom.filename)):
        thiscfg = CFG()
        for thisattr in cfgphantom_attrs:
            setattr(thiscfg.phantom, thisattr, deepcopy(getattr(cfg.phantom, thisattr)[i]))

        cfg_list.append(thiscfg)

    return cfg_list

def PhantomWrapper(cfg):
    if isinstance(cfg.phantom.callback, list):
        # if duplicate phantom callback, raise an error
        if (len(cfg.phantom.callback) != len(set(cfg.phantom.callback))) and len(cfg.phantom.callback) > 1:
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
