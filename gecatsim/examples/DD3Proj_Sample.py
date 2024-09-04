# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import numpy as np
from gecatsim.pyfiles.C_DD3Proj import DD3Proj

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
dzdx = 1
imgXoffset = 0
imgYoffset = 0
imgZoffset = 0
nrviews= 315
viewangles = np.linspace(0, np.pi*2, nrviews, dtype=np.single)  # view angle of each view
zshifts = np.zeros(nrviews, dtype=np.single)  # z-position of each view
nrcols = 51
nrrows = 51
nrplanes = 30
pOrig = np.zeros((nrrows, nrcols, nrplanes), dtype=np.single)
pOrig[25,25,:] = 1
pOrig[25,0,:] = 1
pOrig[0,0,:] = 1

#--------- Run DD3Proj
sinogram = DD3Proj(x0, y0, z0, 
    nrdetcols, nrdetrows, 
    xds, yds, zds, 
    dzdx,
    imgXoffset, imgYoffset, imgZoffset, 
    viewangles, 
    zshifts, 
    nrviews, 
    nrcols, nrrows, nrplanes, 
    pOrig)

#--------- Show results
import matplotlib.pyplot as plt
plt.imshow(sinogram[:,:,7])
plt.show()
