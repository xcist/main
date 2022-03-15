from math import ceil

def mapConfigVariablesToFDK(cfg):

    sid = cfg.scanner.sid
    sdd = cfg.scanner.sdd
    nRow = cfg.scanner.detectorRowsPerMod
    nCol = cfg.scanner.detectorColsPerMod
    nMod = ceil(cfg.scanner.detectorColCount/nCol)
    rowSize = cfg.scanner.detectorRowSize
    colSize = cfg.scanner.detectorColSize
    modWidth = cfg.scanner.detectorColsPerMod*colSize
    dectorYoffset = cfg.scanner.detectorColOffset
    dectorZoffset = cfg.scanner.detectorRowOffset

    fov = cfg.recon.fov
    imageSize = cfg.recon.imageSize
    sliceCount = cfg.recon.sliceCount
    sliceThickness = cfg.recon.sliceThickness
    centerOffset = cfg.recon.centerOffset
    startAngle = cfg.recon.startAngle
    kernelType = cfg.recon.kernelType

    return sid, sdd, nRow, nCol, nMod, rowSize, colSize, modWidth, dectorYoffset, dectorZoffset, \
           fov, imageSize, sliceCount, sliceThickness, centerOffset, startAngle, kernelType