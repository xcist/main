import numpy as np
from skimage.filters import threshold_otsu
from scipy import ndimage

def calculate_mtf_qct(image, fov, image_size, sid=540.0, sdd=950.0, det_pitch=1.0):
    """
    Calculate MTF using the research-based Edge Spread Function (ESF) method.
    Calibrated in lp/mm based on scanner geometry.
    
    Args:
        image: Reconstructed slice (HU or mu)
        fov: Reconstruction field-of-view (mm)
        image_size: Number of pixels (cols)
        sid: Source-to-isocenter distance (mm)
        sdd: Source-to-detector distance (mm)
        det_pitch: Detector element size (mm)
        
    Returns:
        freq: Frequency axis in lp/mm
        mtf: Normalized MTF values
    """
    # 1. Normalize image
    background_roi = image[0:50, 0:50]
    CT_b = np.mean(background_roi)
    CT_o_real = np.max(image)
    image_norm = (image - CT_b) / (CT_o_real - CT_b)
    
    # 2. Find center of the phantom
    thresh = threshold_otsu(image_norm)
    binary = image_norm > thresh
    cy, cx = ndimage.center_of_mass(binary)
    
    # 3. ESF with oversampling
    pixel_size = fov / image_size
    h, w = image.shape
    y, x = np.ogrid[:h, :w]
    r = np.sqrt((x - cx)**2 + (y - cy)**2)
    
    max_r = int(min(cx, cy, w-cx, h-cy))
    r_step = 0.1 # 10x oversampling
    r_bins = np.arange(0, max_r, r_step)
    esf = np.zeros_like(r_bins)
    
    for i, rb in enumerate(r_bins):
        mask = (r >= rb) & (r < rb + r_step)
        if np.any(mask):
            esf[i] = np.mean(image_norm[mask])
        else:
            esf[i] = esf[i-1] if i > 0 else 0.0
    
    # 4. LSF: find the primary edge (focusing on the sharpest transition)
    lsf_full = np.abs(np.diff(esf)) / r_step
    # Skip center to avoid logo artifacts
    peak_idx = np.argmax(lsf_full[int(len(lsf_full)*0.3):]) + int(len(lsf_full)*0.3)
    window = 50 
    lsf_window = lsf_full[max(0, peak_idx-window):min(len(lsf_full), peak_idx+window)]
    
    # Hanning window to reduce FFT leakage
    hann = np.hanning(len(lsf_window))
    lsf_window = lsf_window * hann
    
    # 5. MTF via FFT
    n_fft = 2048
    mtf_complex = np.fft.fft(lsf_window, n=n_fft)
    mtf = np.abs(mtf_complex[:n_fft//2])
    
    if mtf[0] > 0:
        mtf = mtf / mtf[0]
    
    # 6. Frequency axis calibration (lp/mm)
    sample_interval_mm = r_step * pixel_size
    freq = np.arange(n_fft//2) / (n_fft * sample_interval_mm)
    
    # Theoretical limits for information
    # det_nyquist = 1 / (2 * det_pitch * (sid/sdd))
    
    return freq, mtf
