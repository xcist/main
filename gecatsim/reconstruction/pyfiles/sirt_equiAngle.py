"""
SIRT (Simultaneous Iterative Reconstruction Technique) for equiangular geometry
SIRT is similar to CGLS but uses a simpler update rule that tends to be more stable
"""

import numpy as np
from gecatsim.pyfiles.C_DD3Proj import DD3Proj
from gecatsim.pyfiles.C_DD3Back import DD3Back
from gecatsim.reconstruction.pyfiles.mapConfigVariablesToFDK import mapConfigVariablesToFDK

def sirt_equiAngle(cfg, prep, iterations=None, initial_image=None, relaxation=1.0):
    """
    SIRT reconstruction for equiangular CT geometry
    
    Args:
        cfg: CatSim configuration object
        prep: Preprocessed projection data [views, slices, detectors]
        iterations: Number of iterations (default: cfg.recon.num_iterations)
        initial_image: Initial image estimate [rows, cols, slices] in mu space
        relaxation: Relaxation parameter (0 < relaxation <= 1), smaller = more stable
    
    Returns:
        Reconstructed image [rows, cols, slices] in mu space
    """
    
    if iterations is None:
        iterations = getattr(cfg.recon, 'num_iterations', 10)
    
    print("* Using SIRT reconstruction")
    print(f"* Number of iterations: {iterations}")
    print(f"* Relaxation parameter: {relaxation}")
    
    # Get geometry parameters (same as CGLS)
    sid, sdd, nMod, rowSize, modWidth, dectorYoffset, dectorZoffset, \
    fov, imageSize, sliceCount, sliceThickness, centerOffset, startView, rotdir, kernelType \
        = mapConfigVariablesToFDK(cfg)
    
    # Prepare sinogram data
    sino = prep.transpose(0, 2, 1).astype(np.float32)  # [nrviews, nrdetcols, nrdetrows]
    
    nrviews = sino.shape[0]
    nrdetcols = sino.shape[1]
    nrdetrows = sino.shape[2]
    nrcols = imageSize
    nrrows = imageSize
    nrplanes = sliceCount
    
    print(f"* Image size: {nrcols}x{nrrows}x{nrplanes}")
    print(f"* Sinogram size: {nrviews}x{nrdetcols}x{nrdetrows}")
    
    # Source coordinates
    x0 = 0.0
    y0 = sid
    z0 = 0.0
    
    # Detector coordinates
    det_col_indices = np.arange(nrdetcols, dtype=np.float32) - (nrdetcols - 1) * 0.5 + dectorYoffset
    det_row_indices = np.arange(nrdetrows, dtype=np.float32) - (nrdetrows - 1) * 0.5 + dectorZoffset
    xds = det_col_indices * modWidth
    yds = np.full(nrdetcols, sid - sdd, dtype=np.float32)
    zds = det_row_indices * rowSize
    dzdx = 1.0
    
    # View angles
    viewangles = np.linspace(0, 2*np.pi*rotdir, nrviews, endpoint=False, dtype=np.float32)
    
    # Z shifts (for helical)
    zshifts = np.zeros(nrviews, dtype=np.float32)
    
    # Image offsets
    imgXoffset = centerOffset[0]
    imgYoffset = centerOffset[1]
    imgZoffset = centerOffset[2]
    
    # Initialize image
    if initial_image is not None:
        print("* Using provided initial seed")
        img = initial_image.copy().astype(np.float32)
    else:
        print("* Initializing from zero")
        img = np.zeros((nrrows, nrcols, nrplanes), dtype=np.float32)
    
    # Compute normalization factors
    print("* Computing SIRT normalization factors...")
    
    # Column normalization (C): backproject all-ones sinogram
    ones_sino = np.ones_like(sino, dtype=np.float32)
    C = np.zeros((nrrows, nrcols, nrplanes), dtype=np.float32)
    DD3Back(x0, y0, z0,
            nrdetcols, nrdetrows,
            xds, yds, zds,
            dzdx,
            imgXoffset, imgYoffset, imgZoffset,
            viewangles,
            zshifts,
            nrviews,
            ones_sino,
            nrcols, nrrows, nrplanes,
            C)
    C[C < 1e-6] = 1.0  # Avoid division by zero
    
    # Row normalization (R): project all-ones image
    ones_img = np.ones((nrrows, nrcols, nrplanes), dtype=np.float32)
    R = DD3Proj(x0, y0, z0,
                nrdetcols, nrdetrows,
                xds, yds, zds,
                dzdx,
                imgXoffset, imgYoffset, imgZoffset,
                viewangles,
                zshifts,
                nrviews,
                nrcols, nrrows, nrplanes,
                ones_img)
    R[R < 1e-6] = 1.0  # Avoid division by zero
    
    print("* Starting SIRT iterations...")
    
    for i in range(iterations):
        # Forward project current image
        proj = DD3Proj(x0, y0, z0,
                       nrdetcols, nrdetrows,
                       xds, yds, zds,
                       dzdx,
                       imgXoffset, imgYoffset, imgZoffset,
                       viewangles,
                       zshifts,
                       nrviews,
                       nrcols, nrrows, nrplanes,
                       img)
        
        # Compute residual (error)
        residual = sino - proj
        
        # Normalize residual by row sums
        residual_normalized = residual / R
        
        # Backproject normalized residual
        correction = np.zeros((nrrows, nrcols, nrplanes), dtype=np.float32)
        DD3Back(x0, y0, z0,
                nrdetcols, nrdetrows,
                xds, yds, zds,
                dzdx,
                imgXoffset, imgYoffset, imgZoffset,
                viewangles,
                zshifts,
                nrviews,
                residual_normalized,
                nrcols, nrrows, nrplanes,
                correction)
        
        # Normalize correction by column sums and apply relaxation
        correction_normalized = (correction / C) * relaxation
        
        # Update image
        img = img + correction_normalized
        
        # Enforce non-negativity (optional, but common for CT)
        img = np.maximum(img, 0)
        
        # Print progress
        if i == 0 or (i + 1) % 5 == 0 or (i + 1) == iterations:
            residual_norm = np.sqrt(np.mean(residual**2))
            print(f"  Iteration {i+1}: residual = {residual_norm:.2e}")
    
    print("* SIRT reconstruction completed.")
    print(f"* Output: min={np.min(img):.6f}, max={np.max(img):.6f}, mean={np.mean(img):.6f}")
    
    return img
