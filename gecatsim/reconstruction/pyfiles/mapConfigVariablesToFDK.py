from math import ceil
from copy import deepcopy

def mapConfigVariablesToFDK(cfg):

    sid = cfg.scanner.sid
    sdd = cfg.scanner.sdd
    nMod = ceil(cfg.scanner.detectorColCount/cfg.scanner.detectorColsPerMod)
    rowSize = cfg.scanner.detectorRowSize
    modWidth = cfg.scanner.detectorColsPerMod*cfg.scanner.detectorColSize
    dectorYoffset = -cfg.scanner.detectorColOffset
    dectorZoffset = cfg.scanner.detectorRowOffset

    fov = cfg.recon.fov
    imageSize = cfg.recon.imageSize
    sliceCount = cfg.recon.sliceCount
    sliceThickness = cfg.recon.sliceThickness
    centerOffset = deepcopy(cfg.recon.centerOffset)
    # Pass desired X as Y
    centerOffset[1] = deepcopy(cfg.recon.centerOffset[0])
    # Pass desired Y as X; this minus sign is needed to get correct offset
    centerOffset[0] = -deepcopy(cfg.recon.centerOffset[1])

    # The FDK recon seems to be using a "start view" rather than a "start angele".
    # This is a hack until that gets fixed.
    startView_at_view_angle_equals_0 = cfg.protocol.viewCount/2
    if cfg.recon.startAngle <= 180:
        startView = startView_at_view_angle_equals_0 + int(cfg.protocol.viewCount*cfg.recon.startAngle/360)
    elif cfg.recon.startAngle > 180:
        startView = startView_at_view_angle_equals_0 - int(cfg.protocol.viewCount*cfg.recon.startAngle/360)
    else:
        raise Exception("******** Error! Invalid start angle = {} specified. ********".format(cfg.recon.startAngle))

    rotdir = cfg.protocol.rotationDirection
    kernelType = cfg.recon.kernelType

    return sid, sdd, nMod, rowSize, modWidth, dectorYoffset, dectorZoffset, \
           fov, imageSize, sliceCount, sliceThickness, centerOffset, startView, rotdir, kernelType
