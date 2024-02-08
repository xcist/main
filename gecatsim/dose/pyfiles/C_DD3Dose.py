# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/main/blob/master/license/LICENSE

import numpy as np
import os
from ctypes import *
from numpy.ctypeslib import ndpointer
    
def C_DD3Dose(cfg,x0 = None,y0 = None,z0 = None,nrdetcols = None,nrdetrows = None,xds = None,yds = None,zds = None,xoffset = None,yoffset = None,zoffset = None,viewangles = None,zshifts = None,nrviews = None,sino = None,nrcols = None,nrrows = None,nrplanes = None,vol = None,dosevol = None,vox_xy_size = None,vox_z_size = None, debug=True): 
    # was [nrcols, nrrows, nrplanes]
    # convert to [nrplanes, nrrows, nrcols]
    #vol_cp = np.copy(np.transpose(vol, (2,1,0)), order='F')
    vol_cp = np.copy(vol, order='C')
    dosevol_cp = np.copy(dosevol, order='C')
    #dosevol_cp = np.copy(np.transpose(dosevol, (2,1,0)), order='F')
    sinosize = sino.shape
    volsize = dosevol_cp.shape
    xds = np.copy(xds, order='C')
    yds = np.copy(yds, order='C')
    zds = np.copy(zds, order='C')
    #xds = xds[:,None]
    #yds = yds[:,None]
    #zds = zds[None,:]
    #zshifts = zshifts[None,:]
    #viewangles = viewangles[None,:]
    
    if debug:
        ll = cdll.LoadLibrary
        lib_file = "Dose_Recon_Library_Linux64.so"
        recon_lib = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../lib")
        clib = ll(os.path.join(recon_lib, lib_file))
        fun = clib.DD3Dose
    else:
        ###------- load C lib
        cfgnew = cfg.get_current_cfg()
    
        ###------- C function and interface
        fun = cfgnew.dose.doselib.DD3Dose

    fun.argtypes = [c_float, c_float, c_float,
        c_int, c_int,
        ndpointer(c_float), ndpointer(c_float), ndpointer(c_float),
        c_float, c_float, c_float,
        ndpointer(c_float),
        ndpointer(c_float),
        c_int,
        ndpointer(c_float),  # sinogram
        c_int, c_int, c_int,
        ndpointer(c_float),
        ndpointer(c_float),
        c_float, c_float]
    fun.restype = None
    
    ###------- run the function
    sinogram = np.copy(sino.T, order='C')
    #sinogram = np.copy(sino)
   
    fun(x0, y0, z0, 
        nrdetcols, nrdetrows, 
        xds, yds, zds, 
        xoffset, yoffset, zoffset, 
        viewangles, 
        zshifts, 
        nrviews, 
        sinogram, 
        nrcols, nrrows, nrplanes, 
        vol_cp, dosevol_cp, vox_xy_size, vox_z_size);

    dosevol_cp = np.reshape(dosevol_cp, volsize)
    sinogram = np.reshape(sinogram, sinosize)
    # convert to [nrcols, nrrows, nrplanes]
    #dosevol_cp = np.transpose(dosevol_cp,(2,1,0))
    #breakpoint()
    return dosevol_cp, sinogram
    
if __name__ == '__main__':
    from scipy import io as sio
    mat_in = sio.loadmat("dd3dose_in.mat")
    breakpoint()

    pydosevol, pysino = C_DD3Dose(None,mat_in['x0'][0,0], mat_in['y0'][0,0], mat_in['z0'][0,0], mat_in['nrdetcols'][0,0], mat_in['nrdetrows'][0,0], mat_in['xds'], mat_in['yds'], mat_in['zds'],
        mat_in['xoffset'][0,0], mat_in['yoffset'][0,0], mat_in['zoffset'][0,0], mat_in['viewangles'], mat_in['zshifts'], mat_in['nrviews'][0,0], mat_in['sino'], mat_in['nrcols'][0,0], mat_in['nrrows'][0,0], mat_in['nrplanes'][0,0],
        mat_in['vol'], mat_in['dosevol'], mat_in['vox_xy_size'][0,0], mat_in['vox_z_size'][0,0], debug=True)

    mat_out = sio.loadmat("dd3dose_out.mat")

    breakpoint()
    try:
        assert np.allclose(pydosevol, mat_out['dosevol'], atol=1.E-8), 'dosevol'
    except AssertionError as err:
        print(err)
        breakpoint()
