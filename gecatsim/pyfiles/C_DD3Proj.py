# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

from ctypes import *
from numpy.ctypeslib import ndpointer
from gecatsim.pyfiles.CommonTools import *

'''
Wrapper of DD3 Projector.
All vector/matrix input/output are numpy array.
Mingye Wu
'''
def DD3Proj(x0, y0, z0,  # source coordinates
    nrdetcols, nrdetrows,  # detector dimension
    xds, yds, zds,  # detector center coordinates
    dzdx,
    imgXoffset, imgYoffset, imgZoffset, 
    viewangles,  # view angle of each view
    zshifts,  # z-position of each view
    nrviews,  # number of views
    nrcols, nrrows, nrplanes,  # original image 
    pOrig):

    ###------- load C lib
    clib = load_C_lib()
    
    ###------- C function and interface
    fun = clib.DD3Proj
    fun.argtypes = [c_float, c_float, c_float,
        c_int, c_int,
        ndpointer(c_float), ndpointer(c_float), ndpointer(c_float),
        c_float,
        c_float, c_float, c_float,
        ndpointer(c_float),
        ndpointer(c_float),
        c_int,
        ndpointer(c_float),  # sinogram
        c_int, c_int, c_int,
        ndpointer(c_float)]
    fun.restype = None
    
    ###------- run the function
    sinogram = np.zeros([nrviews, nrdetcols, nrdetrows], dtype=np.single)
    
    fun(x0, y0, z0, 
        nrdetcols, nrdetrows, 
        xds, yds, zds, 
        dzdx,
        imgXoffset, imgYoffset, imgZoffset, 
        viewangles, 
        zshifts, 
        nrviews, 
        sinogram, 
        nrcols, nrrows, nrplanes, 
        pOrig)
    
    return sinogram
    
