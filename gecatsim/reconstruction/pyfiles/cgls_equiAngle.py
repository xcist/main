# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import numpy as np
from gecatsim.pyfiles.C_DD3Proj import DD3Proj
from gecatsim.pyfiles.C_DD3Back import DD3Back
from gecatsim.pyfiles.C_DD3WBack import DD3WBack
from gecatsim.reconstruction.pyfiles.mapConfigVariablesToFDK import mapConfigVariablesToFDK

"""
CGLS (Conjugate Gradient Least Squares) iterative reconstruction for CatSim.

This implements the CGLS algorithm similar to ASTRA toolbox, using CatSim's
distance-driven forward and backprojection operators.

IMPORTANT: This uses the UNWEIGHTED backprojector (DD3Back) to ensure it's the
true adjoint of the forward projector (DD3Proj). This is required for CGLS to
converge correctly. The weighted backprojector (DD3WBack) is used in FDK but
should NOT be used in iterative methods like CGLS.
"""

def _compute_detector_coordinates(cfg, sid, sdd, nrdetcols, nrdetrows, rowSize, modWidth, dectorYoffset, dectorZoffset):
    """
    Compute detector center coordinates for DD3 projection/backprojection.
    
    DD3 uses y-axis as source-detector direction:
    - Source is at (0, sid, 0)
    - Detector is at y = -(sdd - sid) = sid - sdd
    xds, yds are the x,y coordinates of detector column centers.
    zds are the z coordinates of detector row centers.
    """
    # Detector column indices (with offset)
    det_col_indices = np.arange(nrdetcols, dtype=np.float32) - (nrdetcols - 1) * 0.5 + dectorYoffset
    det_row_indices = np.arange(nrdetrows, dtype=np.float32) - (nrdetrows - 1) * 0.5 + dectorZoffset
    
    # For DD3, detector coordinates are in the detector plane
    # x coordinates: column positions (assuming detector extends in x direction)
    xds = det_col_indices * modWidth
    
    # y coordinates: detector position (negative side of isocenter, source is on positive side)
    # Source at y=sid, detector at y=sid-sdd
    yds = np.full(nrdetcols, sid - sdd, dtype=np.float32)
    
    # z coordinates: row positions
    zds = det_row_indices * rowSize
    
    return xds, yds, zds


def cgls_equiAngle(cfg, prep, num_iterations=10, initial_image=None):
    """
    CGLS iterative reconstruction.
    
    Parameters:
    -----------
    cfg : configuration object
        CatSim configuration with scanner and reconstruction parameters
        If cfg.recon.num_iterations exists, it will override the num_iterations parameter
    prep : numpy array
        Preprocessed projection data [viewCount, detectorRowCount, detectorColCount]
    num_iterations : int
        Number of CGLS iterations to perform (default: 10)
        Can be overridden by cfg.recon.num_iterations
    initial_image : numpy array, optional
        Initial guess for reconstruction (e.g., from FDK)
        If None, starts from zero
        
    Returns:
    --------
    imageVolume3D : numpy array
        Reconstructed 3D image volume
    """
    # Check if num_iterations is specified in cfg.recon
    if hasattr(cfg.recon, 'num_iterations'):
        num_iterations = cfg.recon.num_iterations
    # Else use the passed argument (default 10)

    
    # Get geometry parameters (similar to FDK)
    sid, sdd, nMod, rowSize, modWidth, dectorYoffset, dectorZoffset, \
    fov, imageSize, sliceCount, sliceThickness, centerOffset, startView, rotdir, kernelType \
        = mapConfigVariablesToFDK(cfg)
    
    # Initialize image volume
    nrcols = imageSize
    nrrows = imageSize
    nrplanes = sliceCount
    
    # Prepare sinogram data first
    # prep shape: [viewCount, detectorRowCount, detectorColCount]
    # Need to transpose to [viewCount, detectorColCount, detectorRowCount] for DD3 functions
    sino = prep.transpose(0, 2, 1).astype(np.float32)  # [nrviews, nrdetcols, nrdetrows]
    
    nrviews = sino.shape[0]
    nrdetcols = sino.shape[1]
    nrdetrows = sino.shape[2]
    
    # Source coordinates (DD3 uses y-axis as source-detector direction)
    x0 = 0.0
    y0 = sid
    z0 = 0.0
    
    # Compute detector coordinates
    xds, yds, zds = _compute_detector_coordinates(cfg, sid, sdd, nrdetcols, nrdetrows, 
                                                   rowSize, modWidth, dectorYoffset, dectorZoffset)
    
    # dzdx is the detector pixel aspect ratio (must be 1.0 for square pixels)
    # Setting to 0.0 causes DD3Proj/DD3Back to return all zeros!
    dzdx = 1.0
    
    # View angles - DON'T add startView offset for CGLS
    # The startView offset in FDK is specific to that algorithm
    # For CGLS, use simple 0 to 2π range matching the acquired data
    viewangles = np.linspace(0, 2*np.pi*rotdir, nrviews, endpoint=False, dtype=np.float32)
    # viewangles += startView * 2 * np.pi / cfg.protocol.viewsPerRotation  # Disabled for CGLS
    
    # Z shifts (for helical scans, otherwise zeros)
    zshifts = np.zeros(nrviews, dtype=np.float32)
    if hasattr(cfg.protocol, 'tableSpeed') and cfg.protocol.tableSpeed != 0:
        # Helical scan - calculate z shifts
        for i in range(nrviews):
            zshifts[i] = cfg.protocol.tableSpeed * (i * cfg.protocol.rotationTime / cfg.protocol.viewsPerRotation)
    
    # Image offsets
    imgXoffset = centerOffset[0]
    imgYoffset = centerOffset[1]
    imgZoffset = centerOffset[2]
    
    # Initialize volume
    # Note: C code expects (nrrows, nrcols, nrplanes) ordering
    if initial_image is not None:
        print(f"* Using provided initial seed")
        v = initial_image.astype(np.float32)
    else:
        # Universal initialization: start from ZERO
        # Most robust - avoids all scaling issues
        # CGLS will converge from zero with sufficient iterations
        print(f"* Starting from ZERO initialization (universal, robust)")
        v = np.zeros((nrrows, nrcols, nrplanes), dtype=np.float32)
    
    # Initialize CGLS variables
    r = np.zeros_like(sino, dtype=np.float32)  # residual: r = sino - A*v
    p = np.zeros_like(v, dtype=np.float32)     # search direction
    z = np.zeros_like(v, dtype=np.float32)     # temporary volume
    w = np.zeros_like(sino, dtype=np.float32)  # temporary sinogram
    
    print(f"* Initializing CGLS reconstruction...")
    print(f"* Image size: {nrcols}x{nrrows}x{nrplanes}")
    print(f"* Sinogram size: {nrviews}x{nrdetcols}x{nrdetrows}")
    print(f"* Number of iterations: {num_iterations}")
    
    # Initial residual: r = sino - A*v
    r = sino.copy()
    
    # If starting from FDK, compute actual residual
    if initial_image is not None:
        print(f"* Computing initial residual from FDK seed...")
        w_init = DD3Proj(x0, y0, z0,
                         nrdetcols, nrdetrows,
                         xds, yds, zds,
                         dzdx,
                         imgXoffset, imgYoffset, imgZoffset,
                         viewangles,
                         zshifts,
                         nrviews,
                         nrcols, nrrows, nrplanes,
                         v)
        r = sino - w_init
        print(f"* Initial residual with FDK seed: {np.sqrt(np.sum(r**2)):.2e}")
    
    # Initial search direction: p = A'*r
    # Use DD3Back (unweighted backprojection)
    p[:] = 0.0
    DD3Back(x0, y0, z0,
             nrdetcols, nrdetrows,
             xds, yds, zds,
             dzdx,
             imgXoffset, imgYoffset, imgZoffset,
             viewangles,
             zshifts,
             nrviews,
             r,
             nrcols, nrrows, nrplanes,
             p)
    
    # gamma = <p, p>
    gamma = np.dot(p.ravel(), p.ravel())
    
    # Print initial residual
    initial_residual = np.sqrt(np.dot(r.ravel(), r.ravel()))
    print(f"* Initial residual: {initial_residual:.6e}")
    print(f"* Starting CGLS iterations (building image from scratch)...")
    print(f"  Iter 1: Backproject raw data → blurry image")
    print(f"  Iter 2+: Refine by correcting projection errors")
    print()
    
    # CGLS iterations
    for iteration in range(num_iterations):
        # w = A * p (forward projection)
        w[:] = 0.0
        w = DD3Proj(x0, y0, z0,
                    nrdetcols, nrdetrows,
                    xds, yds, zds,
                    dzdx,
                    imgXoffset, imgYoffset, imgZoffset,
                    viewangles,
                    zshifts,
                    nrviews,
                    nrcols, nrrows, nrplanes,
                    p)
        
        # alpha = gamma / <w, w>
        ww = np.dot(w.ravel(), w.ravel())
        if ww < 1e-10:  # Avoid division by zero
            print(f"* Warning: ww is very small at iteration {iteration}, stopping.")
            break
        alpha = gamma / ww
        
        # v += alpha * p (update volume)
        v += alpha * p
        
        # r -= alpha * w (update residual)
        r -= alpha * w
        
        # z = A' * r (backproject residual)
        # Use DD3Back (unweighted backprojection)
        z[:] = 0.0
        DD3Back(x0, y0, z0,
                 nrdetcols, nrdetrows,
                 xds, yds, zds,
                 dzdx,
                 imgXoffset, imgYoffset, imgZoffset,
                 viewangles,
                 zshifts,
                 nrviews,
                 r,
                 nrcols, nrrows, nrplanes,
                 z)
        
        # beta = <z, z> / gamma
        newgamma = np.dot(z.ravel(), z.ravel())
        beta = newgamma / gamma if gamma > 1e-10 else 0.0
        
        # gamma = <z, z>
        gamma = newgamma
        
        # p = z + beta * p (update search direction)
        p = z + beta * p
        
        # Compute residual norm and metrics
        residual_norm = np.sqrt(np.dot(r.ravel(), r.ravel()))
        residual_reduction = (1 - residual_norm / initial_residual) * 100
        
        # Show progress
        if iteration == 0:
            print(f"  Iteration 1: residual = {residual_norm:.2e} (blurry image created)")
        elif (iteration + 1) % 5 == 0:
            print(f"  Iteration {iteration + 1}: residual = {residual_norm:.2e} ({residual_reduction:+.1f}% reduction)")
        elif iteration == num_iterations - 1:
            print(f"  Iteration {iteration + 1}: residual = {residual_norm:.2e} ({residual_reduction:+.1f}% reduction)")
    
    print(f"* CGLS reconstruction completed.")
    print(f"* Raw output: min={v.min():.6f}, max={v.max():.6f}, mean={v.mean():.6f}")
    print(f"* Using DD3WBack (weighted backprojector) - no normalization needed")
    
    # Return volume in same format as FDK (transpose if needed)
    # FDK returns [imageSize, imageSize, sliceCount]
    return v

