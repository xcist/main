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
            setattr(thiscfg.phantom, thisattr, deepcopy(getattr(cfg.phantom, thisattr)[i]))

        cfg_list.append(thiscfg)

    return cfg_list


def PhantomWrapper(cfg):
    """
    Prepare phantom configuration for the current view (initialization phase).
    
    ==========================================================================
    WHEN THIS IS CALLED
    ==========================================================================
    Called by RunModels() during the FIRST view/subview only (initialization),
    via the 'Phantom' callback in callback_map. The recalc_map typically has:
        'Phantom': ['viewId', 'subViewId']
    meaning it only recalculates when view or subview changes.
    
    ==========================================================================
    PURPOSE
    ==========================================================================
    Handles DYNAMIC PHANTOM configuration where each view may have different
    phantom data (useful for motion simulation, contrast agent wash-in, etc.).
    
    Detection: Checks if projectorCallback is a NESTED LIST:
        - Nested list (dynamic): [["proj_view0"], ["proj_view1"], ...]
        - Simple list (multiple phantoms): ["proj1", "proj2", ...]
        - String (single phantom): "C_Projector_Voxelized"
    
    ==========================================================================
    WHAT IT DOES FOR DYNAMIC PHANTOMS
    ==========================================================================
    1. Stores original phantom config in cfg._original_dynamic_phantom
       (needed to restore full config for subsequent views)
    
    2. Extracts phantom attributes at index [viewId] from all list attributes
       Example: cfg.phantom.filename[viewId] → temp_cfg.phantom.filename
    
    3. Copies all other cfg attributes to temp_cfg (detector, source, etc.)
    
    4. Returns temp_cfg which now has single-view phantom configuration
    
    ==========================================================================
    FOR NON-DYNAMIC PHANTOMS
    ==========================================================================
    Simply returns cfg unchanged - no processing needed.
    
    Args:
        cfg: Main configuration object with cfg.viewId set
    
    Returns:
        cfg (possibly modified temp_cfg for dynamic phantoms)
    
    Note: The actual phantom data loading (Phantom_Voxelized) and projection
    happen later in ProjectorWrapper, not here. This only prepares the config.
    """
    # ========================================================================
    # DYNAMIC PHANTOM DETECTION
    # ========================================================================
    # Check if projectorCallback is a nested list (list of lists)
    # This indicates different phantom configurations per view
    # Example: [["C_Projector"], ["C_Projector", "Analytic"], ...]
    #          ↑ view 0          ↑ view 1 (multiple phantoms)
    if (isinstance(cfg.phantom.projectorCallback, list) and 
        len(cfg.phantom.projectorCallback) > 0 and 
        isinstance(cfg.phantom.projectorCallback[0], list)):
        
        # ====================================================================
        # PRESERVE ORIGINAL CONFIG
        # ====================================================================
        # Store the full dynamic phantom config so ProjectorWrapper can
        # restore it for subsequent views (each view needs its own config)
        cfg._original_dynamic_phantom = deepcopy(cfg.phantom)
        
        # ====================================================================
        # EXTRACT VIEW-SPECIFIC PHANTOM ATTRIBUTES
        # ====================================================================
        # Get all phantom attributes (filename, projectorCallback, materialList, etc.)
        cfgphantom_attrs = [x for x in dir(cfg.phantom) if not x.startswith('__')]
        
        # Create a temporary cfg to hold the extracted single-view config
        temp_cfg = CFG()
        
        # For each phantom attribute, extract the value at index [viewId]
        # if it's a list with enough elements, otherwise keep as-is
        for thisattr in cfgphantom_attrs:
            phantom_attr = getattr(cfg.phantom, thisattr)
            if isinstance(phantom_attr, list) and len(phantom_attr) > cfg.viewId:
                # Extract view-specific value: e.g., filename[viewId]
                setattr(temp_cfg.phantom, thisattr, phantom_attr[cfg.viewId])
            else:
                # Keep non-list or short-list attributes unchanged
                setattr(temp_cfg.phantom, thisattr, phantom_attr)
        
        # ====================================================================
        # COPY NON-PHANTOM CFG ATTRIBUTES
        # ====================================================================
        # Preserve all other configuration (detector, source, spectrum, etc.)
        for attr in dir(cfg):
            if not attr.startswith('__') and attr != 'phantom':
                setattr(temp_cfg, attr, getattr(cfg, attr))
        
        # Return the view-specific configuration
        cfg = temp_cfg

    # For non-dynamic phantoms, return cfg unchanged
    return cfg

def ProjectorWrapper(cfg, viewId, subViewId):
    """
    Execute phantom projection for the current view and subview.
    
    ==========================================================================
    WHEN THIS IS CALLED
    ==========================================================================
    Called by RunModels() for EVERY view/subview during phantom_scan, via the
    'Projector' callback in callback_map. The recalc_map typically has:
        'Projector': ['viewId', 'subViewId']
    
    This is called AFTER all other RunModels callbacks (detector, source,
    gantry, spectrum, filter, flux) have updated cfg for the current geometry.
    
    ==========================================================================
    INPUT STATE
    ==========================================================================
    At entry, cfg contains:
        - cfg.thisSubView: [totalNumCells, nEbin] float32 array
          Already populated with detection flux (from Mammo_Detection_Flux)
          This represents the X-ray intensity reaching each detector cell
          before passing through the phantom.
    
    ==========================================================================
    OUTPUT STATE
    ==========================================================================
    At exit, cfg.thisSubView has been MULTIPLIED by phantom transmission:
        cfg.thisSubView *= transmission
    
    Where transmission = exp(-sum of path integrals through all materials)
    
    ==========================================================================
    THREE MODES OF OPERATION
    ==========================================================================
    
    MODE 1: DYNAMIC PHANTOM (nested list projectorCallback)
    --------------------------------------------------------
    projectorCallback = [["proj_view0"], ["proj_view1"], ...]
    
    - Restores original phantom config from _original_dynamic_phantom
    - Extracts view-specific phantom attributes at index [viewId]
    - Then continues with Mode 2 or Mode 3 logic
    
    MODE 2: MULTIPLE PHANTOMS PER VIEW (simple list projectorCallback)
    -------------------------------------------------------------------
    projectorCallback = ["C_Projector_Voxelized", "Analytic_Projector", ...]
    
    - Allocates partial_subview: [nSamples, totalNumCells, nEbin]
      ⚠️ WARNING: This can be ~450 GB for full mammo resolution!
    
    - For CPU projector: pre-multiplies by source weights
    - Iterates over each phantom:
        1. SplitCfgPhantom → get individual phantom configs
        2. CopyCfgPhantom → update cfg with this phantom
        3. feval(phantom.callback) → load phantom data (e.g., Phantom_Voxelized)
        4. feval(projectorCallback) → run projector (e.g., C_Projector_Voxelized)
           The projector accumulates into partial_subview
    
    - Combines results:
        CPU: thisSubView *= sum(partial_subview, axis=0)  # sum over samples
        GPU: thisSubView *= partial_subview[0,:,:]        # GPU handles samples internally
    
    - Restores original phantom config for next view
    
    MODE 3: SINGLE PHANTOM (string projectorCallback)
    -------------------------------------------------
    projectorCallback = "C_Projector_Voxelized"
    
    - Simplest case, no partial_subview allocation needed
    - feval(phantom.callback) → load phantom data
    - feval(projectorCallback) → run projector
    - Projector directly modifies cfg.thisSubView
    
    ==========================================================================
    PROJECTOR CALLBACK DETAILS (C_Projector_Voxelized)
    ==========================================================================
    The projector callback performs:
    
    1. Double loop: for srcId in nSamples, for matId in nMaterials
    2. Calls C library: clib.voxelized_projector(phantom, rayOrigins, rayDirections)
    3. Computes path integrals through voxelized phantom
    4. Accumulates transmission: trans += weight * exp(-path_integral)
    5. Final: thisSubView *= trans (or partial_subview for multiple phantoms)
    
    Memory in projector (per call, float64):
        - matPVS: [totalNumCells, nEbin] ~ 100 GB
        - trans: [totalNumCells, nEbin] ~ 100 GB  
        - pValueSpectrum: [totalNumCells, nEbin] ~ 100 GB
    These are freed when the projector function returns.
    
    ==========================================================================
    MEMORY IMPACT
    ==========================================================================
    Mode 2 (multiple phantoms) is the most memory-intensive due to partial_subview.
    Peak memory during multiple-phantom projection:
        partial_subview (450 GB) + projector buffers (300 GB) = ~750 GB
    
    Mode 3 (single phantom) avoids partial_subview:
        Only projector buffers (~300 GB)
    
    Args:
        cfg: Main configuration with thisSubView already containing detection flux
        viewId: Current view index (0 to protocol.viewCount-1)
        subViewId: Current subview index (0 to subViewCount-1)
    
    Returns:
        cfg with cfg.thisSubView multiplied by phantom transmission
    """
    # ==========================================================================
    # MODE 1: DYNAMIC PHANTOM HANDLING
    # ==========================================================================
    # Check if projectorCallback is a nested list (different config per view)
    # This check runs EVERY view to extract the correct phantom for this viewId
    if (isinstance(cfg.phantom.projectorCallback, list) and 
        len(cfg.phantom.projectorCallback) > 0 and 
        isinstance(cfg.phantom.projectorCallback[0], list)):
        
        # Restore the full dynamic phantom config (may have been modified by previous view)
        if hasattr(cfg, '_original_dynamic_phantom'):
            cfg.phantom = deepcopy(cfg._original_dynamic_phantom)
        
        # Extract phantom attributes at index [viewId]
        cfgphantom_attrs = [x for x in dir(cfg.phantom) if not x.startswith('__')]
        
        temp_cfg = CFG()
        
        for thisattr in cfgphantom_attrs:
            phantom_attr = getattr(cfg.phantom, thisattr)
            if isinstance(phantom_attr, list) and len(phantom_attr) > viewId:
                # Extract: filename[viewId], projectorCallback[viewId], etc.
                setattr(temp_cfg.phantom, thisattr, phantom_attr[viewId])
            else:
                setattr(temp_cfg.phantom, thisattr, phantom_attr)
        
        # Preserve non-phantom cfg attributes
        for attr in dir(cfg):
            if not attr.startswith('__') and attr != 'phantom':
                setattr(temp_cfg, attr, getattr(cfg, attr))
        
        # Continue with view-specific phantom config
        cfg = temp_cfg
    
    # ==========================================================================
    # MODE 2: MULTIPLE PHANTOMS PER VIEW
    # ==========================================================================
    # Check if projectorCallback is a simple list (not nested) = multiple phantoms
    if isinstance(cfg.phantom.projectorCallback, list):
        # Split the multi-phantom config into individual phantom configs
        cfg_list = SplitCfgPhantom(cfg)
        
        # ======================================================================
        # ALLOCATE PARTIAL_SUBVIEW - LARGEST MEMORY ALLOCATION IN PIPELINE
        # ======================================================================
        # Shape: [nSamples, totalNumCells, nEbin]
        # This stores transmission for each source sample separately because
        # multiple phantoms may have different spatial relationships to each sample
        # ⚠️ WARNING: ~450 GB for full mammo (9 samples × 170M cells × 74 bins × 4B)
        cfg.partial_subview = np.ones([cfg.src.nSamples, cfg.det.totalNumCells, cfg.spec.nEbin], dtype=np.single)
        
        if not cfg.protocol.gpu_optimization:
            # ================================================================== 
            # CPU PROJECTOR PATH
            # ==================================================================
            # Pre-multiply by source sample weights for proper flux weighting
            # Shape broadcast: [nSamples,1,1] * [nSamples, nCells, nEbin]
            cfg.partial_subview *= cfg.src.weights[:, np.newaxis, np.newaxis]
        
        # ======================================================================
        # ITERATE OVER EACH PHANTOM
        # ======================================================================
        for thiscfg in cfg_list:
            # Copy this phantom's config into main cfg
            cfg = CopyCfgPhantom(thiscfg, cfg)
            
            # Load phantom data into C library (e.g., Phantom_Voxelized)
            # This reads the .raw file and sets up material attenuation tables
            cfg = feval(cfg.phantom.callback, cfg)
            
            # Execute projection (e.g., C_Projector_Voxelized)
            # This multiplies partial_subview by this phantom's transmission
            cfg = feval(cfg.phantom.projectorCallback, cfg, viewId, subViewId)
        
        # ======================================================================
        # COMBINE PARTIAL RESULTS
        # ======================================================================
        # After all phantoms processed, combine into thisSubView
        if not cfg.protocol.gpu_optimization:
            # CPU: Sum over source samples (weighted sum already applied above)
            # Physics: Total transmission = weighted sum of sample transmissions
            cfg.thisSubView *= np.sum(cfg.partial_subview, axis=0)
        else:
            # GPU: GPU projector handles source sample weighting internally
            # Only need the first "sample" which contains combined result
            cfg.thisSubView *= cfg.partial_subview[0,:,:]
        
        # Restore original phantom config for next view (important for dynamic phantom)
        cfg.phantom = deepcopy(cfg._original_dynamic_phantom)
    
    # ==========================================================================
    # MODE 3: SINGLE PHANTOM
    # ==========================================================================
    else:
        # Simplest case: one phantom, one projector
        # No partial_subview needed - projector writes directly to thisSubView
        
        # Load phantom data into C library
        cfg = feval(cfg.phantom.callback, cfg)
        
        # Execute projection - directly modifies cfg.thisSubView
        cfg = feval(cfg.phantom.projectorCallback, cfg, viewId, subViewId)

    return cfg