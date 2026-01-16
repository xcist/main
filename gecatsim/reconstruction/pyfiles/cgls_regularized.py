"""
Regularized CGLS (Conjugate Gradient Least Squares)
Implements Hybrid CGLS with interleaved denoising steps to handle noise retention.
"""

import numpy as np
from scipy.ndimage import gaussian_filter, median_filter
from gecatsim.pyfiles.C_DD3Proj import DD3Proj
from gecatsim.pyfiles.C_DD3Back import DD3Back
from gecatsim.reconstruction.pyfiles.mapConfigVariablesToFDK import mapConfigVariablesToFDK

def _compute_detector_coordinates(cfg, sid, sdd, nrdetcols, nrdetrows, rowSize, modWidth, dectorYoffset, dectorZoffset):
    det_col_indices = np.arange(nrdetcols, dtype=np.float32) - (nrdetcols - 1) * 0.5 + dectorYoffset
    det_row_indices = np.arange(nrdetrows, dtype=np.float32) - (nrdetrows - 1) * 0.5 + dectorZoffset
    xds = det_col_indices * modWidth
    yds = np.full(nrdetcols, sid - sdd, dtype=np.float32)
    zds = det_row_indices * rowSize
    return xds, yds, zds

def cgls_regularized(cfg, prep, num_iterations=10, initial_image=None, 
                     regularization_method=None, regularization_strength=1.0):
    """
    Regularized CGLS reconstruction.
    
    Args:
        cfg: CatSim configuration object
        prep: Preprocessed projection data [views, slices, detectors]
        num_iterations: Number of iterations
        initial_image: Initial guess (mu space)
        regularization_method: 'gaussian', 'median', or None
        regularization_strength: sigma for gaussian, kernel size for median
    """
    
    print(f"* Using Regularized CGLS ({regularization_method}, strength={regularization_strength})")
    
    # Geometry setup
    sid, sdd, nMod, rowSize, modWidth, dectorYoffset, dectorZoffset, \
    fov, imageSize, sliceCount, sliceThickness, centerOffset, startView, rotdir, kernelType \
        = mapConfigVariablesToFDK(cfg)
    
    nrcols, nrrows, nrplanes = imageSize, imageSize, sliceCount
    
    # Sinogram setup [views, cols, rows]
    sino = prep.transpose(0, 2, 1).astype(np.float32)
    nrviews, nrdetcols, nrdetrows = sino.shape
    
    # Coordinates
    x0, y0, z0 = 0.0, sid, 0.0
    xds, yds, zds = _compute_detector_coordinates(cfg, sid, sdd, nrdetcols, nrdetrows, 
                                                   rowSize, modWidth, dectorYoffset, dectorZoffset)
    dzdx = 1.0
    viewangles = np.linspace(0, 2*np.pi*rotdir, nrviews, endpoint=False, dtype=np.float32)
    zshifts = np.zeros(nrviews, dtype=np.float32)
    imgXoffset, imgYoffset, imgZoffset = centerOffset
    
    # Initialize volume
    if initial_image is not None:
        v = initial_image.astype(np.float32)
    else:
        v = np.zeros((nrrows, nrcols, nrplanes), dtype=np.float32)
        
    # Initialize CGLS variables
    r = sino.copy()
    p = np.zeros_like(v)
    z = np.zeros_like(v)
    w = np.zeros_like(sino)
    
    # Initial residual calculation if seeded
    if initial_image is not None:
        print("* Computing initial residual...")
        w_init = DD3Proj(x0, y0, z0, nrdetcols, nrdetrows, xds, yds, zds, dzdx,
                         imgXoffset, imgYoffset, imgZoffset, viewangles, zshifts,
                         nrviews, nrcols, nrrows, nrplanes, v)
        r = sino - w_init
        
    # Initial gradient
    p[:] = 0.0
    DD3Back(x0, y0, z0, nrdetcols, nrdetrows, xds, yds, zds, dzdx,
            imgXoffset, imgYoffset, imgZoffset, viewangles, zshifts,
            nrviews, r, nrcols, nrrows, nrplanes, p)
            
    gamma = np.dot(p.ravel(), p.ravel())
    initial_residual = np.sqrt(np.sum(r**2))
    print(f"* Initial residual: {initial_residual:.2e}")
    
    for i in range(num_iterations):
        # 1. Forward Project
        w[:] = 0.0
        w = DD3Proj(x0, y0, z0, nrdetcols, nrdetrows, xds, yds, zds, dzdx,
                    imgXoffset, imgYoffset, imgZoffset, viewangles, zshifts,
                    nrviews, nrcols, nrrows, nrplanes, p)
                    
        ww = np.sum(w**2)
        if ww < 1e-10: break
        
        alpha = gamma / ww
        
        # 2. Update Image
        v += alpha * p
        
        # 3. Regularization Step (The Fix)
        if regularization_method == 'gaussian':
            # Apply Gaussian smoothing
            # Note: We apply it to the volume v
            # This is a heuristic 'Regularization by Denoising'
            v = gaussian_filter(v, sigma=regularization_strength)
            
        elif regularization_method == 'median':
            v = median_filter(v, size=int(regularization_strength))
            
        # 4. Update Residual
        r -= alpha * w
        
        # 5. Backproject Residual
        z[:] = 0.0
        DD3Back(x0, y0, z0, nrdetcols, nrdetrows, xds, yds, zds, dzdx,
                imgXoffset, imgYoffset, imgZoffset, viewangles, zshifts,
                nrviews, r, nrcols, nrrows, nrplanes, z)
                
        newgamma = np.sum(z**2)
        beta = newgamma / gamma if gamma > 1e-10 else 0.0
        gamma = newgamma
        
        p = z + beta * p
        
        # Progress
        res_norm = np.sqrt(np.sum(r**2))
        if i==0 or (i+1)%5==0 or i==num_iterations-1:
            print(f"  Iter {i+1}: residual={res_norm:.2e}")
            
    return v
