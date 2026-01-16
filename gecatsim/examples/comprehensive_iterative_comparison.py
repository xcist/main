import gecatsim as xc
import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import time
from skimage.filters import threshold_otsu
from scipy import ndimage

print("="*80)
print("Comprehensive Iterative Reconstruction Comparison (Master Suite)")
print("ART vs SIRT vs CGLS - Continuous runs with refined QCT-MTF (lp/mm)")
print("="*80)

# 1. Configuration and Data Loading
# ============================================================================
ct = xc.CatSim("./cfg/Phantom_Sample", "./cfg/Scanner_Sample_generic", "./cfg/Protocol_Sample_axial")
ct.resultsName = "test_compare"

# Simulation Parameters
mu_water = 0.02
hu_offset = -1000
snapshot_iters = [100, 300, 500]
slice_idx = 2

# Distance-Driven Geometry setup
from gecatsim.pyfiles.C_DD3Proj import DD3Proj
from gecatsim.pyfiles.C_DD3Back import DD3Back
from gecatsim.reconstruction.pyfiles.mapConfigVariablesToFDK import mapConfigVariablesToFDK

sid, sdd, nMod, rowSize, modWidth, dectorYoffset, dectorZoffset, \
fov, imageSize, sliceCount, sliceThickness, centerOffset, startView, rotdir, kernelType \
    = mapConfigVariablesToFDK(ct)

# Check for existing data or run
if not os.path.exists("test_compare.prep"):
    print("Generating noisy projection data...")
    ct.run_all()

print("Loading data...")
sino = xc.rawread(ct.resultsName + ".prep",
                  [ct.protocol.viewCount, ct.scanner.detectorRowCount, ct.scanner.detectorColCount],
                  'float')

img_phantom = xc.rawread('test_compare_TRUE_PHANTOM_512x512x4.raw', [4, 512, 512], 'float')
img_fdk = xc.rawread('test_compare_FDK_512x512x4.raw', [4, 512, 512], 'float')

nrviews, nrdetrows, nrdetcols = sino.shape
nrcols = ct.recon.imageSize
nrrows = ct.recon.imageSize
nrplanes = ct.recon.sliceCount

# Geometry variables for DD3
x0, y0, z0 = 0.0, sid, 0.0
det_col_indices = np.arange(nrdetcols) - (nrdetcols - 1) * 0.5 + dectorYoffset
det_row_indices = np.arange(nrdetrows) - (nrdetrows - 1) * 0.5 + dectorZoffset
xds = det_col_indices.astype(np.float32) * modWidth
yds = np.full(nrdetcols, sid - sdd, dtype=np.float32)
zds = det_row_indices.astype(np.float32) * rowSize
dzdx = 1.0
viewangles = np.linspace(0, 2*np.pi*rotdir, nrviews, endpoint=False, dtype=np.float32)
zshifts = np.zeros(nrviews, dtype=np.float32)
imgXoffset, imgYoffset, imgZoffset = centerOffset

results = {}

# 2. Iterative Reconstructions
# ============================================================================

# --- SIRT ---
print("\nRunning SIRT...")
results['SIRT'] = {}
v_sirt = np.zeros((nrrows, nrcols, nrplanes), dtype=np.float32)
ones_sino = np.ones_like(sino); C_sirt = np.zeros((nrrows, nrcols, nrplanes), dtype=np.float32)
DD3Back(x0, y0, z0, nrdetcols, nrdetrows, xds, yds, zds, dzdx, imgXoffset, imgYoffset, imgZoffset, viewangles, zshifts, nrviews, ones_sino, nrcols, nrrows, nrplanes, C_sirt)
C_sirt[C_sirt < 1e-6] = 1.0
ones_img = np.ones((nrrows, nrcols, nrplanes), dtype=np.float32)
R_sirt = DD3Proj(x0, y0, z0, nrdetcols, nrdetrows, xds, yds, zds, dzdx, imgXoffset, imgYoffset, imgZoffset, viewangles, zshifts, nrviews, nrcols, nrrows, nrplanes, ones_img)
R_sirt[R_sirt < 1e-6] = 1.0

for i in range(500):
    proj = DD3Proj(x0, y0, z0, nrdetcols, nrdetrows, xds, yds, zds, dzdx, imgXoffset, imgYoffset, imgZoffset, viewangles, zshifts, nrviews, nrcols, nrrows, nrplanes, v_sirt)
    residual_normalized = (sino - proj) / R_sirt
    correction = np.zeros((nrrows, nrcols, nrplanes), dtype=np.float32)
    DD3Back(x0, y0, z0, nrdetcols, nrdetrows, xds, yds, zds, dzdx, imgXoffset, imgYoffset, imgZoffset, viewangles, zshifts, nrviews, residual_normalized, nrcols, nrrows, nrplanes, correction)
    v_sirt += 0.5 * (correction / C_sirt); v_sirt = np.maximum(v_sirt, 0)
    if (i+1) in snapshot_iters:
        results['SIRT'][i+1] = {'image': (v_sirt.transpose(2,0,1)*(1000/mu_water)+hu_offset).copy()}

# --- ART ---
print("Running ART...")
results['ART'] = {}
v_art = np.zeros((nrrows, nrcols, nrplanes), dtype=np.float32)
subset_size = 10; indices = list(range(0, nrviews, subset_size))
for it in range(500):
    np.random.shuffle(indices)
    for start_idx in indices:
        end_idx = min(start_idx + subset_size, nrviews)
        sub_sino = sino[start_idx:end_idx]; sub_angles = viewangles[start_idx:end_idx]; sub_zshifts = zshifts[start_idx:end_idx]
        sub_proj = DD3Proj(x0, y0, z0, nrdetcols, nrdetrows, xds, yds, zds, dzdx, imgXoffset, imgYoffset, imgZoffset, sub_angles, sub_zshifts, len(sub_angles), nrcols, nrrows, nrplanes, v_art)
        sub_R = DD3Proj(x0, y0, z0, nrdetcols, nrdetrows, xds, yds, zds, dzdx, imgXoffset, imgYoffset, imgZoffset, sub_angles, sub_zshifts, len(sub_angles), nrcols, nrrows, nrplanes, np.ones_like(v_art))
        sub_C = np.zeros((nrrows, nrcols, nrplanes), dtype=np.float32)
        DD3Back(x0, y0, z0, nrdetcols, nrdetrows, xds, yds, zds, dzdx, imgXoffset, imgYoffset, imgZoffset, sub_angles, sub_zshifts, len(sub_angles), np.ones_like(sub_sino), nrcols, nrrows, nrplanes, sub_C)
        correction = np.zeros((nrrows, nrcols, nrplanes), dtype=np.float32)
        DD3Back(x0, y0, z0, nrdetcols, nrdetrows, xds, yds, zds, dzdx, imgXoffset, imgYoffset, imgZoffset, sub_angles, sub_zshifts, len(sub_angles), (sub_sino - sub_proj)/np.maximum(sub_R, 1e-6), nrcols, nrrows, nrplanes, correction)
        v_art += 0.5 * (correction / np.maximum(sub_C, 1e-6)); v_art = np.maximum(v_art, 0)
    if (it+1) in snapshot_iters:
        results['ART'][it+1] = {'image': (v_art.transpose(2,0,1)*(1000/mu_water)+hu_offset).copy()}

# --- CGLS ---
print("Running CGLS...")
results['CGLS'] = {}
v_cgls = np.zeros((nrrows, nrcols, nrplanes), dtype=np.float32); r = sino.copy(); p = np.zeros_like(v_cgls)
DD3Back(x0, y0, z0, nrdetcols, nrdetrows, xds, yds, zds, dzdx, imgXoffset, imgYoffset, imgZoffset, viewangles, zshifts, nrviews, r, nrcols, nrrows, nrplanes, p)
gamma = np.sum(p**2)
for i in range(500):
    w = DD3Proj(x0, y0, z0, nrdetcols, nrdetrows, xds, yds, zds, dzdx, imgXoffset, imgYoffset, imgZoffset, viewangles, zshifts, nrviews, nrcols, nrrows, nrplanes, p)
    ww = np.sum(w**2); alpha = gamma / ww if ww > 1e-10 else 0; v_cgls += alpha * p; r -= alpha * w
    z = np.zeros_like(v_cgls); DD3Back(x0, y0, z0, nrdetcols, nrdetrows, xds, yds, zds, dzdx, imgXoffset, imgYoffset, imgZoffset, viewangles, zshifts, nrviews, r, nrcols, nrrows, nrplanes, z)
    newgamma = np.sum(z**2); p = z + (newgamma / gamma if gamma > 1e-10 else 0) * p; gamma = newgamma
    if (i+1) in snapshot_iters:
        results['CGLS'][i+1] = {'image': (v_cgls.transpose(2,0,1)*(1000/mu_water)+hu_offset).copy()}

# 3. Refined MTF Calculation (QCT Method)
# ============================================================================
def calculate_mtf_qct_paper_method(image, fov, image_size):
    background_roi = image[0:50, 0:50]; CT_b = np.mean(background_roi); image_norm = (image - CT_b) / (np.max(image) - CT_b)
    thresh = threshold_otsu(image_norm); binary = image_norm > thresh; cy, cx = ndimage.center_of_mass(binary)
    pixel_size = fov / image_size; h, w = image.shape; y, x = np.ogrid[:h, :w]; r = np.sqrt((x - cx)**2 + (y - cy)**2)
    max_r = int(min(cx, cy, w-cx, h-cy)); r_step = 0.1; r_bins = np.arange(0, max_r, r_step); esf = np.zeros_like(r_bins)
    for i, rb in enumerate(r_bins):
        mask = (r >= rb) & (r < rb + r_step)
        if np.any(mask): esf[i] = np.mean(image_norm[mask])
        else: esf[i] = esf[i-1] if i > 0 else 0.0
    lsf_full = np.abs(np.diff(esf)) / r_step; peak_idx = np.argmax(lsf_full[int(len(lsf_full)*0.3):]) + int(len(lsf_full)*0.3)
    window = 50; lsf_window = lsf_full[max(0, peak_idx-window):min(len(lsf_full), peak_idx+window)]
    lsf_window *= np.hanning(len(lsf_window)); n_fft = 2048; mtf_complex = np.fft.fft(lsf_window, n=n_fft); mtf = np.abs(mtf_complex[:n_fft//2])
    if mtf[0] > 0: mtf = mtf / mtf[0]
    return np.arange(n_fft//2) / (n_fft * r_step * pixel_size), mtf

# 4. Final Plotting
# ============================================================================
print("\nGenerating final comparison figure...")
phantom_slice = img_phantom[slice_idx, :, :]; fdk_slice = img_fdk[slice_idx, :, :]
freq_axis, phantom_mtf = calculate_mtf_qct_paper_method(phantom_slice, fov, imageSize)
_, fdk_mtf = calculate_mtf_qct_paper_method(fdk_slice, fov, imageSize)

fig = plt.figure(figsize=(30, 20)); gs = GridSpec(5, 5, figure=fig, hspace=0.35, wspace=0.25)
for row_idx, n_iter in enumerate(snapshot_iters):
    if row_idx == 0:
        ax = fig.add_subplot(gs[0, 0]); ax.imshow(phantom_slice, cmap='gray', vmin=-200, vmax=200); ax.set_title('Phantom (GT)', fontweight='bold'); ax.axis('off')
        ax = fig.add_subplot(gs[0, 1]); ax.imshow(fdk_slice, cmap='gray', vmin=-200, vmax=200); ax.set_title(f'FDK (Baseline)\nRMSE: {np.sqrt(np.mean((fdk_slice-phantom_slice)**2)):.1f} HU', fontweight='bold'); ax.axis('off')
    else:
        ax = fig.add_subplot(gs[row_idx, 0:2]); ax.text(0.5, 0.5, f'{n_iter} Iterations', ha='center', va='center', fontsize=20, fontweight='bold', transform=ax.transAxes); ax.axis('off')
    for col_idx, method in enumerate(['SIRT', 'ART', 'CGLS']):
        ax = fig.add_subplot(gs[row_idx, col_idx+2]); img = results[method][n_iter]['image'][slice_idx, :, :]
        ax.imshow(img, cmap='gray', vmin=-200, vmax=200); ax.set_title(f'{method} ({n_iter} iter)\nRMSE: {np.sqrt(np.mean((img-phantom_slice)**2)):.1f} HU'); ax.axis('off')

ax_mtf = fig.add_subplot(gs[3, :]); ax_mtf.plot(freq_axis, phantom_mtf, 'k--', linewidth=3, label='Phantom')
ax_mtf.plot(freq_axis, fdk_mtf, 'gray', linewidth=2.5, label='FDK', alpha=0.7)
colors = {'ART':'red','SIRT':'blue','CGLS':'green'}; linestyles = {100:'-', 300:'--', 500:':'}
for m in ['SIRT', 'ART', 'CGLS']:
    for ni in snapshot_iters:
        _, m_val = calculate_mtf_qct_paper_method(results[m][ni]['image'][slice_idx,:,:], fov, imageSize)
        ax_mtf.plot(freq_axis, m_val, color=colors[m], linestyle=linestyles[ni], linewidth=1.5, label=f'{m} {ni}', alpha=0.7)
ax_mtf.set_xlim(0, 1.0); ax_mtf.set_ylim(0, 1.1); ax_mtf.set_xlabel('Spatial Frequency (lp/mm)'); ax_mtf.set_ylabel('MTF'); ax_mtf.legend(ncol=5); ax_mtf.grid(True, alpha=0.3)
ax_mtf.axvline(x=1/(2*1.0*(540/950)), color='r', linestyle=':', alpha=0.5, label='Det Nyquist'); ax_mtf.axvline(x=0.512, color='b', linestyle=':', alpha=0.5, label='Recon Nyquist')

ax_h = fig.add_subplot(gs[4, 0:2]); ax_h.plot(phantom_slice[256, :], 'k--', linewidth=2.5); ax_h.plot(fdk_slice[256, :], 'gray')
ax_v = fig.add_subplot(gs[4, 3:5]); ax_v.plot(phantom_slice[:, 256], 'k--', linewidth=2.5); ax_v.plot(fdk_slice[:, 256], 'gray')
for m in ['SIRT', 'ART', 'CGLS']:
    for ni in snapshot_iters:
        ax_h.plot(results[m][ni]['image'][slice_idx, 256, :], color=colors[m], linestyle=linestyles[ni], alpha=0.5)
        ax_v.plot(results[m][ni]['image'][slice_idx, :, 256], color=colors[m], linestyle=linestyles[ni], alpha=0.5)
ax_h.set_title('Horizontal Profile'); ax_h.grid(True, alpha=0.3); ax_v.set_title('Vertical Profile'); ax_v.grid(True, alpha=0.3)

plt.savefig('iterative_comparison_refined_scaled.png', dpi=150, bbox_inches='tight')
print(f"\n✓ Master Script Complete. Saved figure to iterative_comparison_refined_scaled.png")
