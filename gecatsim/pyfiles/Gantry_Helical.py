# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import math
import copy
from gecatsim.pyfiles.CommonTools import *

def Gantry_Helical(cfg, viewId):
    ### calculate transform matrices
    srcTransform, srcRotation, detTransform, detRotation = helical_transform(cfg, viewId)
    
    ### apply transform
    src = cfg.src
    det = cfg.det
        
    cfg.srcNew = copy.copy(src)
    cfg.srcNew.corners = np.c_[src.corners, np.ones((src.nCorners, 1), dtype=np.single)] @ srcTransform
    cfg.srcNew.samples = np.c_[src.samples, np.ones((src.nSamples, 1), dtype=np.single)] @ srcTransform
    cfg.srcNew.front = src.front @ srcRotation
    cfg.srcNew.lateral = src.lateral @ srcRotation
    cfg.srcNew.long = src.long @ srcRotation
    
    cfg.detNew = copy.copy(det)
    cfg.detNew.modCoords = np.c_[det.modCoords, np.ones((det.nMod, 1), dtype=np.single)] @ detTransform
    cfg.detNew.uvecs = det.uvecs @ detRotation
    cfg.detNew.vvecs = det.vvecs @ detRotation
    
    return cfg

def helical_transform(cfg, viewId):
    # wobble and the 4x4 matrix
    wobble = (viewId%2*2-1)*cfg.protocol.wobbleDistance  # odd and even
    Rwobble = np.array([[1, 0, 0, wobble], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    
    # rotation angle and the 4x4 matrix
    beta = cfg.protocol.startAngle*math.pi/180+cfg.time/cfg.protocol.rotationTime*math.tau*cfg.protocol.rotationDirection
    Rbeta = np.array([[np.cos(beta), -np.sin(beta), 0, 0], [np.sin(beta), np.cos(beta), 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    
    # tilt angle and the 4x4 matrix
    phi = cfg.protocol.tiltAngle*math.pi/180
    Rphi = np.array([[1, 0, 0, 0], [0, np.cos(phi), np.sin(phi), 0], [0, -np.sin(phi), np.cos(phi), 0], [0, 0, 0, 1]])
    
    # translation distance and the 4x4 matrix
    translation = -(cfg.protocol.startZ + cfg.time*cfg.protocol.tableSpeed)
    Rtrans = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, translation], [0, 0, 0, 1]])
    
    # The src transform is Z*phi*beta*wobble
    srcTransform = np.single(Rtrans @ Rphi @ Rbeta @ Rwobble)
    srcTransform = srcTransform[0:3, :].T
    srcRotation = srcTransform[0:3, :]
    
    # The detector transform is Z*phi*beta
    detTransform = np.single(Rtrans @ Rphi @ Rbeta)
    detTransform = detTransform[0:3, :].T
    detRotation = detTransform[0:3, :]
    
    return srcTransform, srcRotation, detTransform, detRotation
