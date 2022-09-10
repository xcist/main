# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

import numpy.matlib as nm
import math
from catsim.pyfiles.CommonTools import *

def Source_Uniform(cfg):
    '''
    Returns a source structure (src) with a uniform focal spot on a plane with a given
    target angle and with a specified width and optical length.
    Mingye Wu, GE Research
    
    '''
    
    # shortcuts
    nx = cfg.physics.srcXSampleCount
    ny = cfg.physics.srcYSampleCount
    nz = ny
    dx = cfg.scanner.focalspotWidth/nx
    dy = -cfg.scanner.focalspotLength/np.tan(cfg.scanner.targetAngle/180*math.pi)/ny
    dz = cfg.scanner.focalspotLength/nz
    
    # sample coords and weights
    nSamples = nx*ny
    weights = nm.repmat(1/nSamples, 1, nSamples)
    
    x = (np.arange(0, nx)-(nx-1)/2)*dx
    y = (np.arange(0, ny)-(ny-1)/2)*dy+cfg.scanner.sid
    z = (np.arange(0, nz)-(nz-1)/2)*dz
    
    x = nm.repmat(x, 1, ny).T
    y = nm.repmat(y, nx, 1).T.reshape(nSamples, 1)
    z = nm.repmat(z, nx, 1).T.reshape(nSamples, 1)
    
    samples = np.c_[x, y, z]
    
    # focal spot offset
    if cfg.protocol.focalspotOffset:
        samples = samples + nm.repmat(cfg.protocol.focalspotOffset, nSamples, 1)
    
    # corners
    if nx==1 and ny==1:
        corners = samples
    elif nx>1 and ny>1:
        corners = np.c_[samples[0, :], samples[nx-1, :], samples[(ny-1)*nx, :], samples[-1, :]].T
    else:
        corners = np.c_[samples[0, :], samples[-1, :]].T
    nCorners = corners.shape[0]
    
    # source definition
    if not cfg.src:
        cfg.src = CFG()
    cfg.src.nSamples = nSamples
    cfg.src.samples = np.single(samples)
    cfg.src.weights = np.single(weights)
    cfg.src.front   = np.array([[0, -1, 0]], dtype=np.single)
    cfg.src.lateral = np.array([[1, 0, 0]], dtype=np.single)
    cfg.src.long    = np.array([[0, 0, 1]], dtype=np.single)
    cfg.src.nCorners = nCorners
    cfg.src.corners = np.single(corners)
    
    return cfg


if __name__ == "__main__":

    cfg = source_cfg("./cfg/default.cfg")
    
    '''
    cfg.scanner.sid = 541
    cfg.scanner.focalspotWidth = 1
    cfg.scanner.focalspotLength = 0.9
    cfg.scanner.targetAngle = 7
    cfg.protocol.focalspotOffset = [0, 0, 0];
    cfg.physics.srcXSampleCount = 4
    cfg.physics.srcYSampleCount = 3
    '''
    
    cfg = Source_Uniform(cfg)
    check_value(cfg.src.samples)
    check_value(cfg.src.corners)
    check_value(cfg.src.lateral)
