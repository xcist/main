"""
ART (Algebraic Reconstruction Technique) / Simultaneous ART (SART) for equiangular geometry
This implementation supports Ordered Subsets (OS-SART) for faster convergence.
"""

import numpy as np
from gecatsim.pyfiles.C_DD3Proj import DD3Proj
from gecatsim.pyfiles.C_DD3Back import DD3Back
from gecatsim.reconstruction.pyfiles.mapConfigVariablesToFDK import mapConfigVariablesToFDK

def art_equiAngle(cfg, prep, iterations=None, initial_image=None, relaxation=0.5, subset_size=10):
    """
    ART/SART reconstruction for equiangular CT geometry with Ordered Subsets.
    
    Args:
        cfg: CatSim configuration object (uses cfg.recon for defaults)
        prep: Preprocessed projection data [views, slices, detectors]
        iterations: Number of full iterations (default: cfg.recon.num_iterations)
        initial_image: Initial image estimate [rows, cols, slices] in mu space
        relaxation: Relaxation parameter (default: 0.5)
        subset_size: Number of views per subset/incremental update (default: 10)
    
    Returns:
        Reconstructed image [rows, cols, slices] in mu space
    """
    
    if iterations is None:
        iterations = getattr(cfg.recon, 'num_iterations', 10)
    
    print("* Using ART (OS-SART) reconstruction")
    print(f"* Total full iterations: {iterations}")
    print(f"* Subset size (views per update): {subset_size}")
    print(f"* Relaxation: {relaxation}")
    
    # Get geometry parameters
    sid, sdd, nMod, rowSize, modWidth, dectorYoffset, dectorZoffset, \
    fov, imageSize, sliceCount, sliceThickness, centerOffset, startView, rotdir, kernelType \
        = mapConfigVariablesToFDK(cfg)
    
    # Prepare sinogram data
    # DD3 expects [views, dets, slices]
    sino = prep.transpose(0, 2, 1).astype(np.float32)  # [nrviews, nrdetcols, nrdetrows]
    
    nrviews = sino.shape[0]
    nrdetcols = sino.shape[1]
    nrdetrows = sino.shape[2]
    nrcols = imageSize
    nrrows = imageSize
    nrplanes = sliceCount
    
    # Coordinate Setup (DD3 style)
    x0, y0, z0 = 0.0, sid, 0.0
    det_col_indices = np.arange(nrdetcols, dtype=np.float32) - (nrdetcols - 1) * 0.5 + dectorYoffset
    det_row_indices = np.arange(nrdetrows, dtype=np.float32) - (nrdetrows - 1) * 0.5 + dectorZoffset
    xds = det_col_indices * modWidth
    yds = np.full(nrdetcols, sid - sdd, dtype=np.float32)
    zds = det_row_indices * rowSize
    dzdx = 1.0
    
    viewangles = np.linspace(0, 2*np.pi*rotdir, nrviews, endpoint=False, dtype=np.float32)
    zshifts = np.zeros(nrviews, dtype=np.float32)
    
    imgXoffset, imgYoffset, imgZoffset = centerOffset
    
    # Initialization
    if initial_image is not None:
        img = initial_image.copy().astype(np.float32)
    else:
        img = np.zeros((nrrows, nrcols, nrplanes), dtype=np.float32)
    
    # Setup subsets
    indices = list(range(0, nrviews, subset_size))
    
    for it in range(iterations):
        # Shuffle subsets for better performance in OS-SART
        np.random.shuffle(indices)
        
        for start_idx in indices:
            end_idx = min(start_idx + subset_size, nrviews)
            sub_sino = sino[start_idx:end_idx, :, :]
            sub_angles = viewangles[start_idx:end_idx]
            sub_zshifts = zshifts[start_idx:end_idx]
            current_subset_views = end_idx - start_idx
            
            # 1. Project current image for this subset
            sub_proj = DD3Proj(x0, y0, z0, nrdetcols, nrdetrows, xds, yds, zds, dzdx,
                               imgXoffset, imgYoffset, imgZoffset, sub_angles, sub_zshifts,
                               current_subset_views, nrcols, nrrows, nrplanes, img)
            
            # 2. Compute residual
            sub_residual = sub_sino - sub_proj
            
            # 3. Compute row-wise normalization (Weighted SART style)
            # R = A * 1_volume
            ones_vol = np.ones((nrrows, nrcols, nrplanes), dtype=np.float32)
            sub_R = DD3Proj(x0, y0, z0, nrdetcols, nrdetrows, xds, yds, zds, dzdx,
                            imgXoffset, imgYoffset, imgZoffset, sub_angles, sub_zshifts,
                            current_subset_views, nrcols, nrrows, nrplanes, ones_vol)
            sub_R[sub_R < 1e-6] = 1.0
            sub_residual_norm = sub_residual / sub_R
            
            # 4. Compute column-wise normalization
            # C = A^T * 1_sino
            ones_sub_sino = np.ones_like(sub_sino)
            sub_C = np.zeros((nrrows, nrcols, nrplanes), dtype=np.float32)
            DD3Back(x0, y0, z0, nrdetcols, nrdetrows, xds, yds, zds, dzdx,
                    imgXoffset, imgYoffset, imgZoffset, sub_angles, sub_zshifts,
                    current_subset_views, ones_sub_sino, nrcols, nrrows, nrplanes, sub_C)
            sub_C[sub_C < 1e-6] = 1.0
            
            # 5. Backproject normalized residual
            correction = np.zeros((nrrows, nrcols, nrplanes), dtype=np.float32)
            DD3Back(x0, y0, z0, nrdetcols, nrdetrows, xds, yds, zds, dzdx,
                    imgXoffset, imgYoffset, imgZoffset, sub_angles, sub_zshifts,
                    current_subset_views, sub_residual_norm, nrcols, nrrows, nrplanes, correction)
            
            # 6. Apply update with relaxation
            img += relaxation * (correction / sub_C)
            
            # Non-negativity constraint
            img = np.maximum(img, 0)
            
        # Logging progress
        if (it + 1) % 5 == 0 or (it + 1) == iterations:
            # Full projection error check (optional, but keep it quiet)
            print(f"  Completed iteration {it + 1}/{iterations}")

    print("* ART reconstruction completed.")
    return img
