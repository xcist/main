# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import sys
import numpy.matlib as nm
from scipy import interpolate, io
import matplotlib.pyplot as plt
from gecatsim.pyfiles.CommonTools import *
from collections import defaultdict, OrderedDict
from mpl_toolkits.axes_grid1 import make_axes_locatable

def GetDefaultWidthLength(shape):
    shape = shape.lower()
    if shape=="uniform":
        width = 1.
        length = 1.
    elif shape=="gaussian":
        width = 1.
        length = 1.
    elif shape=="performix":
        width = 1.
        length = 1.
        #scanner.fs_performix_width = 0.92
        #scanner.fs_performix_length = 0.76
    elif shape=="pharos_small":
        width = 1.
        length = 1.
    elif shape=="pharos_large":
        width = 1.
        length = 1.
    elif shape=="gemini_small":
        width = 1.
        length = 1.
    elif shape=="gemini_large":
        width = 1.
        length = 1.

    return width, length

def GetIntensity(cfg):
    if cfg.physics.srcXSampleCount%2==0:
        nx = 10*cfg.physics.srcXSampleCount
    else:
        nx = 10*cfg.physics.srcXSampleCount+1
    if cfg.physics.srcYSampleCount%2==0:
        ny = 10*cfg.physics.srcYSampleCount
    else:
        ny = 10*cfg.physics.srcYSampleCount+1
    x_grid, y_grid = np.mgrid[0:nx, 0:ny]
    if cfg.scanner.focalspotShape.lower() == "gaussian":
        # one sigma in units of pixel number
        if hasattr(cfg.scanner, "focalspotSigmaX"): sx = cfg.scanner.focalspotSigmaX
        else: sx = 10
        if hasattr(cfg.scanner, "focalspotSigmaZ"): sz = cfg.scanner.focalspotSigmaZ
        else: sz = 10
        weights = np.exp(-((x_grid-(nx-1)/2)**2/sx**2+(y_grid-(ny-1)/2)**2/sz**2)/2)
    elif cfg.scanner.focalspotShape.lower() == 'uniform':
        weights = np.ones((nx, ny))
    weights /= np.sum(weights)

    return weights

def ParseFocalspotData(path):
    _, ext = os.path.splitext(path)
    # ending in .mat
    if '.mat' == ext:
        alldata = io.loadmat(path)
        data = alldata['I'].T
        pixsize_x = alldata['dx'][0,0]
        pixsize_z = alldata['dz'][0,0]
    # ending in .npz
    elif '.npz' == ext:
        alldata = np.load(path)
        data = alldata['data']
        pixsize_x = alldata['pixsize_x']
        pixsize_z = alldata['pixsize_z']

    return data, pixsize_x, pixsize_z, 0, 0

def SetFocalspot(cfg):
    # if shape and data is not defined, defaults to Uniform; 
    if (not hasattr(cfg.scanner, "focalspotShape")) and (not hasattr(cfg.scanner, "focalspotData")):
        cfg.scanner.focalspotShape = "Uniform"
    elif hasattr(cfg.scanner, "focalspotShape") and hasattr(cfg.scanner, "focalspotData"):
    #elif cfg.scanner.focalspotShape and hasattr(cfg.scanner, "focalspotData"):
        print("Focal spot: FocalspotData is set and will override FocalspotShape.")
        delattr(cfg.scanner, "focalspotShape")
        #sys.exit()
    # load default width and length
    if not all([hasattr(cfg.scanner, "focalspotWidth"), hasattr(cfg.scanner, "focalspotLength")]):
        cfg.scanner.focalspotWidth, cfg.scanner.focalspotLength = GetDefaultWidthLength(cfg.scanner.focalspotShape)

    # load npz focus spot image, the measured intensity will always be in the xz plane
    if hasattr(cfg.scanner, "focalspotData"):
        cfg.scanner.focalspotData = my_path.find("focal_spot", cfg.scanner.focalspotData, '')
        I, pixsize_x, pixsize_z, xstart, zstart = ParseFocalspotData(cfg.scanner.focalspotData)
        cfg.scanner.focalspotPixSizeX = pixsize_x
        cfg.scanner.focalspotPixSizeZ = pixsize_z
    else:
        I = GetIntensity(cfg)
        xstart, zstart = 0, 0

    nx, nz = I.shape
    ny = nz
    if hasattr(cfg.scanner, 'focalspotShape'):# and cfg.scanner.focalspotShape.lower() == 'uniform':
        dx = cfg.scanner.focalspotWidth/nx
        dy = -cfg.scanner.focalspotLength/np.tan(cfg.scanner.targetAngle*np.pi/180.)/ny
        dz = cfg.scanner.focalspotLength/nz
    else:
        dx = pixsize_x
        dy = -pixsize_z/np.tan(cfg.scanner.targetAngle*np.pi/180)
        dz = pixsize_z

    # remove too small values
    I /= np.max(I)
    fs_pos_x = (xstart + dx*(np.arange(nx)+0.5))
    fs_pos_z = (zstart + dz*(np.arange(nz)+0.5))
    nx, nz = I.shape

    # recenter, based on com
    valid_ctr_idx = np.where(I>=0.1) # this is only for focal spot center, only calculate com based on intensity >= 10%max
    _idx_x_min = np.min(valid_ctr_idx[0])
    _idx_x_max = np.max(valid_ctr_idx[0])+1
    _idx_z_min = np.min(valid_ctr_idx[1])
    _idx_z_max = np.max(valid_ctr_idx[1])+1

    fs_pos_x -= np.average(fs_pos_x[_idx_x_min:_idx_x_max], weights=np.sum(I[_idx_x_min:_idx_x_max, _idx_z_min:_idx_z_max], axis=1))
    fs_pos_z -= np.average(fs_pos_z[_idx_z_min:_idx_z_max], weights=np.sum(I[_idx_x_min:_idx_x_max, _idx_z_min:_idx_z_max], axis=0))
    
    # rescale
    if not hasattr(cfg.scanner, 'focalspotShape') or cfg.scanner.focalspotShape.lower() != 'uniform':
    #if not cfg.scanner.focalspotShape or cfg.scanner.focalspotShape.lower() != 'uniform':
        Ix = np.sum(I, axis=1) # xis along the vertical axis (axis0), so we need to sum along axis 1
        Ix /= np.max(Ix)
        max_idx = np.argmax(Ix)
        pos1 = np.interp(cfg.scanner.focalspotWidthThreshold, Ix[0:max_idx], fs_pos_x[0:max_idx])
        pos2 = np.interp(cfg.scanner.focalspotWidthThreshold, Ix[max_idx:][::-1], fs_pos_x[max_idx:][::-1])
        W0 = np.abs(pos2 - pos1) # in units of mm

        Iz = np.sum(I, axis=0)
        Iz /= np.max(Iz)
        max_idx = np.argmax(Iz)
        zpos1 = np.interp(cfg.scanner.focalspotLengthThreshold, Iz[0:max_idx], fs_pos_z[0:max_idx])
        zpos2 = np.interp(cfg.scanner.focalspotLengthThreshold, Iz[max_idx:][::-1], fs_pos_z[max_idx:][::-1])
        L0 = np.abs(zpos2 - zpos1)

    # down sampling to match oversampling
    os_nx = cfg.physics.srcXSampleCount
    os_ny = cfg.physics.srcYSampleCount
    os_nz = os_ny
    def GetRange(pos, profile, th):
        '''
        clever way of determining sampling range based on position, profile, and threshold
        '''
        profile /= np.max(profile)
        idx = np.where(profile>th)[0][0]
        if idx==0:
            new_begin = pos[0]
        else:
            new_begin = np.interp(th, profile[idx-1:idx+1], pos[idx-1:idx+1])
        idx = np.where(profile>th)[0][-1]
        if idx==len(pos)-1:
            new_end = pos[-1]
        else:
            new_end = np.interp(th, profile[idx:idx+2][::-1], pos[idx:idx+2][::-1])

        return new_begin, new_end

    # clever way of sampling
    if not hasattr(cfg.scanner, 'focalspotShape') or cfg.scanner.focalspotShape.lower() != 'uniform':
    #if not cfg.scanner.focalspotShape or cfg.scanner.focalspotShape.lower() != 'uniform':
        os_range_x = GetRange(fs_pos_x, np.sum(I, axis=1), 0.02)
        os_range_z = GetRange(fs_pos_z, np.sum(I, axis=0), 0.02)
        os_dx = (os_range_x[1] - os_range_x[0])/os_nx
        os_dz = (os_range_z[1] - os_range_z[0])/os_nz
    else:
        # need to make sure that the os range is consistent with Source_Uniform
        os_range_x = [0.5*(-cfg.scanner.focalspotWidth), None]
        os_range_z = [0.5*(-cfg.scanner.focalspotLength), None]
        os_dx = cfg.scanner.focalspotWidth/os_nx
        os_dz = cfg.scanner.focalspotLength/os_nz
    os_x = os_range_x[0] + (np.arange(os_nx)+0.5)*os_dx # should use start+offset because it may be not symmetric in the range
    os_z = os_range_z[0] + (np.arange(os_nz)+0.5)*os_dz # 0.5 comes from the actual center is not at edge
    os_y = -os_z/np.tan(cfg.scanner.targetAngle*np.pi/180.)
    [os_xx, os_zz] = np.meshgrid(os_x, os_z)
    #os_interp = interpolate.interp2d(fs_pos_z, fs_pos_x, I, kind='linear')

    #os_I = os_interp(os_z, os_x)
    os_interp = interpolate.RectBivariateSpline(fs_pos_z, fs_pos_x, I.T, kx=1, ky=1)
    os_I = os_interp(os_z, os_x).T

    if hasattr(cfg.scanner, 'focalspotData') or cfg.scanner.focalspotShape.lower() != 'uniform':
        os_xx *= cfg.scanner.focalspotWidth/W0
        os_zz *= cfg.scanner.focalspotLength/L0

    os_yy = -os_zz/np.tan(cfg.scanner.targetAngle*np.pi/180.);

    # remove low-weight sampling and normalize
    weights = os_I.T.flatten()
    nSamples = weights.size
    samples = np.c_[os_xx.flatten(), cfg.scanner.sid+os_yy.flatten(), os_zz.flatten()]
    weights /= np.sum(weights)

    # re-center samples based on center of mass
    samples[:,0] -= np.average(samples[:,0], weights=weights)
    samples[:,1] -= np.average(samples[:,1], weights=weights)
    samples[:,1] += cfg.scanner.sid
    samples[:,2] -= np.average(samples[:,2], weights=weights)

    # offset
    if hasattr(cfg.protocol, 'focalspotOffset'):
        samples = samples + nm.repmat(cfg.protocol.focalspotOffset, nSamples, 1)

    # find corners
    if os_nx==1 and os_ny==1:
        corners = samples
    elif os_nx>1 and os_ny>1:
        corners = np.array([[os_x[0], cfg.scanner.sid+os_y[0], os_z[0]],
                            [os_x[0], cfg.scanner.sid+os_y[-1], os_z[-1]],
                            [os_x[-1], cfg.scanner.sid+os_y[-1], os_z[-1]],
                            [os_x[-1], cfg.scanner.sid+os_y[0], os_z[0]]])
    else:
        corners = np.c_[samples[0, :], samples[-1, :]].T
    nCorners = corners.shape[0]
    
    # source definition
    #if not hasattr(cfg, 'src'):
    if not cfg.src:
        cfg.src = CFG()
    cfg.src.nSamples = nSamples
    cfg.src.samples = np.single(samples)
    cfg.src.weights = np.single(weights[None])
    cfg.src.front   = np.array([[0, -1, 0]], dtype=np.single)
    cfg.src.lateral = np.array([[1, 0, 0]], dtype=np.single)
    cfg.src.long    = np.array([[0, 0, 1]], dtype=np.single)
    cfg.src.nCorners = nCorners
    cfg.src.corners = np.single(corners)

    return cfg
