gecatsim\pyfiles\PhantomProjectorWrapper.py
@@ -1,466 +1,94 @@
"""
PhantomProjectorWrapper.py - Phantom and Projector Configuration Management

=============================================================================
OVERVIEW - WHERE THIS FITS IN THE XCIST PIPELINE
=============================================================================

This module is called during the phantom_scan phase of CatSim.run_all():

    CatSim.run_all()
      └── phantom_scan()
            └── one_scan()
                  └── RunModels() with callback_map including:
                        - 'Phantom': PhantomWrapper()     [initialization phase]
                        - 'Projector': ProjectorWrapper() [per-view projection]

The module handles THREE phantom configuration scenarios:

1. SINGLE PHANTOM (simplest case):
   cfg.phantom.projectorCallback = "C_Projector_Voxelized"
   → One phantom, same for all views

2. MULTIPLE PHANTOMS PER VIEW (e.g., breast + lesion + calcifications):
   cfg.phantom.projectorCallback = ["C_Projector_Voxelized", "Analytic_Projector", ...]
   cfg.phantom.filename = ["breast.json", "lesion.json", ...]
   → Multiple phantoms combined, same set for all views

3. DYNAMIC PHANTOM (nested list - different phantom per view):
   cfg.phantom.projectorCallback = [
       ["C_Projector_Voxelized"],           # view 0
       ["C_Projector_Voxelized"],           # view 1
       ["Analytic_Projector", "C_Proj..."], # view 2 (multiple phantoms)
       ...
   ]
   → Each view can have different phantom configuration (e.g., motion simulation)

=============================================================================
MEMORY CONSIDERATIONS
=============================================================================

For multiple phantoms per view, partial_subview is allocated:
    Shape: [nSamples, totalNumCells, nEbin]
    Size: nSamples × totalNumCells × nEbin × 4 bytes
    Example: 9 × 170,572,500 × 74 × 4 = ~450 GB for full mammo resolution!

This is the LARGEST memory allocation in the entire pipeline and occurs
because each phantom's transmission must be computed independently before
combining (multiplicative combination for overlapping phantoms).

For single phantom, no partial_subview is needed - projection writes
directly to cfg.thisSubView.

=============================================================================
FUNCTION CALL SEQUENCE
=============================================================================

Per-view execution flow:

    RunModels(cfg, viewId, subViewId)
        │
        ├── [First view only] PhantomWrapper(cfg)
        │       └── Handles dynamic phantom: extracts cfg.phantom for viewId
        │
        └── ProjectorWrapper(cfg, viewId, subViewId)
                │
                ├── [Dynamic phantom] Re-extract cfg.phantom for viewId
                │
                ├── [Multiple phantoms] 
                │   ├── SplitCfgPhantom() → list of individual phantom cfgs
                │   ├── Allocate partial_subview [nSamples, nCells, nEbin]
                │   ├── Loop over each phantom:
                │   │   ├── CopyCfgPhantom() → copy phantom config
                │   │   ├── feval(phantom.callback) → Phantom_Voxelized()
                │   │   └── feval(projectorCallback) → C_Projector_Voxelized()
                │   └── Combine: thisSubView *= sum(partial_subview)
                │
                └── [Single phantom]
                    ├── feval(phantom.callback) → Phantom_Voxelized()
                    └── feval(projectorCallback) → C_Projector_Voxelized()

=============================================================================
"""

import sys
from copy import copy, deepcopy
from gecatsim.pyfiles.CommonTools import *

def CopyCfgPhantom(cfgfrom, cfgto=None):
    """
    Copy phantom configuration from one CFG object to another.
    
    Used during multiple-phantom iteration to update the main cfg with
    each individual phantom's configuration while preserving all other
    cfg attributes (detector, source, spectrum, etc.).
    
    Args:
        cfgfrom: Source CFG object containing phantom to copy
        cfgto: Destination CFG object (created if None)
    
    Returns:
        cfgto with phantom attributes deep-copied from cfgfrom
    
    Note: Uses deepcopy to avoid reference issues when phantom attributes
    are modified during projection.
    """
    if cfgto is None:
        cfgto = CFG()
    cfgto.phantom = deepcopy(cfgfrom.phantom)

    return cfgto

def SplitCfgPhantom(cfg):
    """
    Split a multi-phantom configuration into a list of single-phantom CFGs.
    
    When cfg.phantom contains lists (multiple phantoms to combine), this
    function creates separate CFG objects for each phantom so they can be
    processed individually by their respective projectors.
    
    Example input:
        cfg.phantom.filename = ["breast.raw", "lesion.raw", "calc.raw"]
        cfg.phantom.projectorCallback = ["C_Projector", "C_Projector", "Analytic"]
        cfg.phantom.materialList = [["adipose"], ["tumor"], ["calcium"]]
    
    Output:
        cfg_list[0].phantom.filename = "breast.raw"
        cfg_list[0].phantom.projectorCallback = "C_Projector"
        cfg_list[0].phantom.materialList = ["adipose"]
        ...and so on for each phantom
    
    Args:
        cfg: Main CFG with phantom attributes as lists
    
    Returns:
        List of CFG objects, each with single phantom configuration
    
    Note: Only phantom attributes are copied; other cfg attributes (det, src, etc.)
    are accessed from the main cfg via CopyCfgPhantom during iteration.
    """
    cfg_list = []
    cfgphantom_attrs = [x for x in dir(cfg.phantom) if not x.startswith('__')]
    for i in range(len(cfg.phantom.filename)):
        thiscfg = CFG()
        for thisattr in cfgphantom_attrs:
            # broadcast in this case
            if not isinstance(getattr(cfg.phantom, thisattr), list) or len(getattr(cfg.phantom, thisattr))==1:
                setattr(thiscfg.phantom, thisattr, deepcopy(getattr(cfg.phantom, thisattr)))
            else:
                setattr(thiscfg.phantom, thisattr, deepcopy(getattr(cfg.phantom, thisattr)[i]))

        cfg_list.append(thiscfg)

    return cfg_list

def PhantomWrapper(cfg):
    # nested loop, like [[pro1, proj2], [proj1], [proj2], ...]
    is_dynamic_phantom = False
    if isinstance(cfg.phantom.callback, list) and isinstance(cfg.phantom.callback[0], list):
        is_dynamic_phantom = True
        _cfg_list = SplitCfgPhantom(cfg)
        _origCfgPhantom = CopyCfgPhantom(cfg)
        cfg = CopyCfgPhantom(_cfg_list[cfg.viewId], cfg)

    if isinstance(cfg.phantom.callback, list):
        # if duplicate phantom callback, raise an error
        if (len(cfg.phantom.callback) != len(set(cfg.phantom.callback))) and len(cfg.phantom.callback) > 1:
            print("Error! XCIST does not support repeated phantom types.\n")
            sys.exit(1)

        cfg_list = SplitCfgPhantom(cfg)
        origCfgPhantom = CopyCfgPhantom(cfg)
        cfg_nom_list = []
        for thiscfg in cfg_list:
            # use cfg here because we want to use cfg's other attr than phantom without the heavy copy overhead
            cfg = CopyCfgPhantom(thiscfg, cfg)
            cfg = feval(cfg.phantom.callback, cfg)
            if hasattr(cfg.phantom, 'numberOfMaterials'):
                cfg_nom_list.append(cfg.phantom.numberOfMaterials)
                delattr(cfg.phantom, 'numberOfMaterials')
            else:
                # ncat does not have such property
                cfg_nom_list.append(None)

        # restore original CFG
        if is_dynamic_phantom:
            cfg = CopyCfgPhantom(_origCfgPhantom, cfg)
            cfg.phantom.numberOfMaterials = [None]*len(cfg.phantom.callback)
            cfg.phantom.numberOfMaterials[cfg.viewId] = cfg_nom_list
        else:
            cfg = CopyCfgPhantom(origCfgPhantom, cfg)
            cfg.phantom.numberOfMaterials = cfg_nom_list
    else:
        cfg = feval(cfg.phantom.callback, cfg)


    return cfg

def ProjectorWrapper(cfg, viewId, subViewId):
    is_dynamic_phantom = False
    if isinstance(cfg.phantom.callback, list) and isinstance(cfg.phantom.callback[0], list):
        is_dynamic_phantom = True
        _cfg_list = SplitCfgPhantom(cfg)
        _origCfgPhantom = CopyCfgPhantom(cfg)
        cfg = CopyCfgPhantom(_cfg_list[cfg.viewId], cfg)

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

    if is_dynamic_phantom:
        cfg = CopyCfgPhantom(_origCfgPhantom, cfg)

    return cfg
