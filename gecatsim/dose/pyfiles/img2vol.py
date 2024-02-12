# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/main/blob/master/license/LICENSE

import os
import numpy as np
from scipy import interpolate
import gecatsim as xc
from gecatsim.dose.pyfiles.xyfovimg import xyfovimg
    
def img2vol(cfg = None): 
    ## read recon image
    if hasattr(cfg.dose, 'imageFileName'):
        if os.path.exists(cfg.dose.imageFileName):
            filename = cfg.dose.imageFileName
        else:
            filename = cfg.dose.imageFileName+'.recon'
    else:
        filename = cfg.results_basename+'.recon'
    
    img = xc.rawread(filename, [cfg.recon.sliceCount, cfg.recon.imageSize, cfg.recon.imageSize], 'float')
    img = np.transpose(img, (2,1,0))
    #img = np.reshape(img, tuple(np.array([cfg.recon_size,cfg.recon_size,cfg.recon_planes])), order="F")
    ## rebin img
    n_voxel = cfg.dose.nVoxel
    img_rebin = np.zeros((n_voxel,n_voxel,cfg.recon.sliceCount))
    for ii in np.arange(cfg.recon.sliceCount):
        X =np.arange(cfg.recon.imageSize)
        Y = np.arange(cfg.recon.imageSize)
        Z = img[:,:,ii]
        XI = np.arange(0,cfg.recon.imageSize, cfg.recon.imageSize / n_voxel)
        YI=  np.arange(0,cfg.recon.imageSize, cfg.recon.imageSize / n_voxel)
        f_interp = interpolate.interp2d(X,Y,Z,kind='linear')
        ZI = f_interp(XI, YI)
        img_rebin[:,:,ii] = ZI

    vol = np.copy(img_rebin)
    ## convert to Mu, in 1/cm
    if cfg.recon.unit.lower() == 'hu' and np.max(vol) > 10:
        vol = (vol - cfg.recon.huOffset) / 1000 * cfg.recon.mu*10 # in cm 
    
    ## remove outer air
    outer_ind = xyfovimg(n_voxel,n_voxel,cfg.recon.sliceCount,n_voxel/2,n_voxel/2,n_voxel/2-0.5,n_voxel/2-0.5)
    outer_ind = outer_ind == 0
    vol[np.logical_or(vol<cfg.dose.outerAirThreshold, outer_ind)] = 0
    # # adjust bone Mu, because of BHC
    # vol(vol>cfg.waterThreshold)=vol(vol>cfg.waterThreshold)*1.04;
    
    tmpmask = np.logical_or(vol<0, np.isnan(vol), np.isinf(vol))
    vol[tmpmask] = 0.0
    vol = np.single(vol)
    ## Density
    mask_bone = vol > cfg.dose.waterThreshold
    mask_soft = 1 - mask_bone
    dens_vol = vol / (mask_soft * cfg.dose.muWater / 1.0 + mask_bone * cfg.dose.muBone / 1.92)
    
    mass_vol = dens_vol * (cfg.recon.fov / 10 / cfg.dose.nVoxel) ** 3

    return vol,dens_vol,mass_vol
