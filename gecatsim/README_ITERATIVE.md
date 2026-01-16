# CatSim Iterative Reconstruction & Analysis Suite

This documentation summarizes the iterative reconstruction algorithms and analysis tools added to the `gecatsim` package.

## 1. Iterative Reconstruction Algorithms

The following algorithms are located in `gecatsim/reconstruction/pyfiles/` and can be called by setting `cfg.recon.reconType` to their respective names.

### Algebraic Reconstruction Technique (ART / OS-SART)
- **Module**: `art_equiAngle.py`
- **Algorithm**: Ordered Subset SART (Simultaneous Algebraic Reconstruction Technique).
- **Key Features**: 
    - Fast convergence using subsets of projection data.
    - Robust handling of noisy raw projections.
    - Integrated non-negativity constraints.
- **Usage**:
  ```python
  ct.recon.reconType = "art_equiAngle"
  ct.recon.num_iterations = 100
  ct.recon.relaxation = 0.5
  img = ct.recon_direct()
  ```

### Simultaneous Iterative Reconstruction Technique (SIRT)
- **Module**: `sirt_equiAngle.py`
- **Algorithm**: Standard SIRT with column/row sum normalization.
- **Key Features**: 
    - Extremely stable convergence.
    - Excellent denoising performance from raw data.
- **Usage**:
  ```python
  ct.recon.reconType = "sirt_equiAngle"
  ct.recon.num_iterations = 300
  img = ct.recon_direct()
  ```

### Conjugate Gradient Least Squares (CGLS)
- **Module**: `cgls_equiAngle.py`
- **Algorithm**: Minimizes $\|Ax - b\|^2$ using conjugate gradients.
- **Key Features**: 
    - Fast mathematical convergence.
    - Strictly minimizes the least-squares error.
- **Usage**:
  ```python
  ct.recon.reconType = "cgls_equiAngle"
  ct.recon.num_iterations = 50
  img = ct.recon_direct()
  ```

---

## 2. Advanced Resolution Analysis (MTF)

The package now includes research-based MTF analysis tools in `gecatsim/pyfiles/mtf_tools.py`.

### QCT Edge-based MTF
- **Method**: Calculates the Modulation Transfer Function (MTF) using the Edge Spread Function (ESF) $\to$ Line Spread Function (LSF) pipeline.
- **Calibration**: Automatically converts frequency units to **lp/mm** by considering:
    - Reconstruction FOV and pixel size.
    - Scanner magnification ($SID/SDD$).
    - Detector pitch.
- **Robustness**: Uses Otsu thresholding for center-finding and Focuses on sharp edges (e.g., circular boundaries) to maintain accuracy even in non-radially-symmetric phantoms.

#### Example Usage:
```python
from gecatsim.pyfiles.mtf_tools import calculate_mtf_qct

freq, mtf = calculate_mtf_qct(image_slice, fov=500, image_size=512)
plt.plot(freq, mtf)
plt.xlabel("Spatial Frequency (lp/mm)")
```

---

## 3. Comprehensive Comparison Script
The script `gecatsim/examples/comprehensive_iterative_comparison.py` provides a turnkey solution for comparing these methods against a TRUE phantom ground truth, including visual comparisons, RMSE metrics, and profile analysis.
