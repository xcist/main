# Copyright 2025, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import numpy as np
from gecatsim.pyfiles.C_DD3Proj_mm import DD3Proj_mm

# assumptions
# voxel size is 0.5 in x/y direction, 2.0 in z direction
# detector size is 1.0 mm x 1.0 mm
vox_xy_size = np.single(0.5)
vox_z_size = np.single(2)

#--------- Set parameters
# source coordinates
x0 = 0
y0 = 550
z0 = 0

# detector dimension and center coordinates
nrdetcols = 200
nrdetrows = 16
xds = np.linspace(-nrdetcols/2+0.5, nrdetcols/2-0.5, nrdetcols, dtype=np.single)
yds = -400*np.ones(nrdetcols, dtype=np.single)
zds = np.linspace(-nrdetrows/2+0.5, nrdetrows/2-0.5, nrdetrows, dtype=np.single)

# original image and view settings
imgXoffset = 0
imgYoffset = 0
imgZoffset = 0
nrviews= 315
viewangles = np.linspace(0, np.pi*2, nrviews, dtype=np.single)  # view angle of each view
zshifts = np.zeros(nrviews, dtype=np.single)  # z-position of each view
nrcols = 51
nrrows = 51
nrplanes = 30
pOrig = np.zeros((nrrows, nrcols, nrplanes), dtype=np.single)  # voxelized object, where each voxel's value is the volume fraction
pOrig[25,25,:] = 1
pOrig[25,0,:] = 1
pOrig[0,0,:] = 1

# a 2D mask to determine if a voxel is included in projection. 1 means included
xy_mask = np.ones((2*nrrows, nrcols), dtype=np.uint8)

#--------- Run DD3Proj
sinogram = DD3Proj_mm(x0, y0, z0, 
    nrdetcols, nrdetrows, 
    xds, yds, zds, 
    imgXoffset, imgYoffset, imgZoffset, 
    viewangles, 
    zshifts, 
    nrviews, 
    nrcols, nrrows, nrplanes, 
    pOrig,
    vox_xy_size, vox_z_size,
    xy_mask
    )

#--------- Show results
import matplotlib.pyplot as plt
plt.imshow(sinogram[:,:,7])
plt.show()
