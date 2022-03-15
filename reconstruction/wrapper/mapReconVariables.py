from catsim.CommonTools import *


def mapReconVariables(cfg):
    fov = cfg.recon.fov
    imagesize = cfg.recon.imageSize
    slicecount = cfg.recon.sliceCount
    kerneltype = cfg.recon.kernelType
    sliceThickness = cfg.recon.sliceThickness
    centerOffset = cfg.recon.centerOffset
    startangle = cfg.recon.startangle

    return fov, imagesize, slicecount, kerneltype, sliceThickness, centerOffset, startangle