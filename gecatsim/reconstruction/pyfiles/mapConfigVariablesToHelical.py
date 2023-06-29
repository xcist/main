from math import ceil
from copy import deepcopy


def mapConfigVariablesToHelical(cfg):

    
    sid = cfg.scanner.sid # source to isocenter distance, in mm
    sdd = cfg.scanner.sdd # source to detector distance, in mm
    YL = int(cfg.scanner.detectorColCount)
    ZL = int(cfg.scanner.detectorRowCount) #
    N_2pi =  cfg.protocol.viewsPerRotation #
    ViewN = cfg.protocol.viewCount
    N_Turn = (cfg.protocol.viewCount-1)/cfg.protocol.viewsPerRotation
    h = cfg.protocol.tableSpeed*cfg.protocol.rotationTime # in xcist, in unit of mm/s
    startAngle = cfg.protocol.startAngle
    dectorYoffset = -cfg.scanner.detectorColOffset   # coloffset is in unit of cols, could be float
    dectorZoffset = cfg.scanner.detectorRowOffset

    # The following lines are used to define the reconstruction paramters
    k1 = 5  # The order to define the 3D weighting function
    delta = 60  # The range to define smoothness of 2D weigthing function
    HSCoef = 0.6  # This is used to define the half-scan range

    # nMod = ceil(cfg.scanner.detectorColCount/cfg.scanner.detectorColsPerMod)
    rowSize = cfg.scanner.detectorRowSize  # size of 1 detector row in unit of mm, float
    ColSize = cfg.scanner.detectorColSize  # size of 1 detector col in unit of mm, float

    imageSize = cfg.recon.imageSize # how many pixels, currently a single int, has to be square
    sliceCount = cfg.recon.sliceCount # how many slices in recon for 3d recon
    sliceThickness = cfg.recon.sliceThickness # float in mm
    objR =0.5* cfg.recon.fov #FOV in mm

    kernelType = cfg.recon.kernelType
    centerOffset = deepcopy(cfg.recon.centerOffset)
    # Pass desired X as Y
    centerOffset[1] = deepcopy(cfg.recon.centerOffset[0])
    # Pass desired Y as X
    centerOffset[0] = -deepcopy(cfg.recon.centerOffset[1])
    # kernelType = cfg.recon.kernelType
    # centerOffset = deepcopy(cfg.recon.centerOffset)
    # # Pass desired X as Y
    # centerOffset[1] = deepcopy(cfg.recon.centerOffset[0])
    # # Pass desired Y as -X
    # centerOffset[0] = -deepcopy(cfg.recon.centerOffset[1])

    return  sid, sdd, YL, ZL, ViewN, N_Turn, N_2pi, h, startAngle, dectorYoffset, dectorZoffset, \
            k1, delta, HSCoef, rowSize, ColSize, imageSize,sliceCount, sliceThickness, centerOffset, objR, kernelType


