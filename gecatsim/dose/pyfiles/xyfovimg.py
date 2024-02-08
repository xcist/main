# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/main/blob/master/license/LICENSE

import numpy as np
import numpy.matlib
    
def xyfovimg(nrrows = None,nrcols = None,nrplanes = 1,radius1 = None,radius2 = None,centerrow = None,centercol = None): 
    
    # Set default input arguments if not provided
    if radius1 is None:
       radius1 = nrrows / 2.0
    
    if radius2 is None:
       radius2 = nrcols / 2.0
    
    if centerrow is None:
       centerrow = (-1 + nrrows) / 2.0
    
    if centercol is None :
        centercol = (-1 + nrcols) / 2.0
    
    # Calculate squared row and column indices ; then find where radius is smaller than input radius
    row2img = (np.arange(nrrows)[:,None] - centerrow) ** 2 * np.ones((1,nrcols))
    col2img = (np.arange(nrcols)[None,:] - centercol) ** 2 * np.ones((nrrows,1))
    img = (row2img / (radius1 ** 2) + col2img / (radius2 ** 2)) <= 1
    
    # Replicate nrplanes times if 3D volume (cylinder) is requested
    if nrplanes > 1:
        #img = np.matlib.repmat(img, np.array([1,1,nrplanes]))
        img = np.repeat(img[:, :, np.newaxis], nrplanes, axis=2)
    return img
