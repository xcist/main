# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

from catsim.CommonTools import *
import matplotlib.pyplot as plt
# Need to import new recons as they are added
from reconstruction.pyfiles.fdk_equiAngle import fdk_equiAngle


def recon(cfg):

    prep = load_prep(cfg)

    # The following line doesn't work - need to fix it when new recons are added.
    # imageVolume3D = feval("reconstruction." + cfg.recon.reconType, cfg, prep)

    # A hack until the previous line is fixed.
    imageVolume3D = fdk_equiAngle(cfg, prep)
    imageVolume3D = scaleReconData(cfg, imageVolume3D)

    if cfg.recon.saveImageVolume:
        saveImageVolume(cfg, imageVolume3D)

    if cfg.recon.saveSingleImages:
        saveSingleImages(cfg, imageVolume3D)

    if cfg.recon.displayImagePictures:
        cfg = displayImagePictures(cfg, imageVolume3D)

    if cfg.recon.saveImagePictureFiles:
        cfg = saveImagePictureFiles(cfg, imageVolume3D)

    return cfg

     
def load_prep(cfg):

    print("* Loading the projection data...")
    prep = rawread(cfg.resultsName + ".prep",
                  [cfg.protocol.viewCount, cfg.scanner.detectorRowCount, cfg.scanner.detectorColCount],
                  'float')
                  
    return prep


def scaleReconData(cfg, imageVolume3D):

    print('* Scaling recon data...')
    if cfg.recon.unit =='HU':
        imageVolume3D = imageVolume3D*(1000/(cfg.recon.mu))
        imageVolume3D = imageVolume3D + cfg.recon.huOffset
    elif cfg.recon.unit == '/mm':
        imageVolume3D = imageVolume3D
    elif cfg.recon.unit == '/cm':
        imageVolume3D = imageVolume3D*10
    else:
        raise Exception('******** Error! An unsupported recon unit was specified: {:s}. ********'.format(cfg.recon.unit))


    return imageVolume3D

    
def saveImageVolume(cfg, imageVolume3D):

    print('* Writing the recon results to one big file...')

    imageVolume3D_size_string = str(cfg.recon.imageSize) + 'x' + str(cfg.recon.imageSize) + 'x' + str(cfg.recon.sliceCount)
    fname = cfg.resultsName + '_' + imageVolume3D_size_string + '.raw'
    imageVolume3D = imageVolume3D.transpose(2, 0, 1)
    imageVolume3D = imageVolume3D.copy(order='C')
    rawwrite(fname, imageVolume3D)


def saveSingleImages(cfg, imageVolume3D):

    print('* Writing the recon results to individual files...')

    sliceIndicesToSave = range(0, cfg.recon.sliceCount)
    for sliceIndexToSave in sliceIndicesToSave:
        imageVolume3D_size_string = str(cfg.recon.imageSize) + 'x' + str(cfg.recon.imageSize) + 'x1'
        sliceNumberString = 'slice' + str(sliceIndexToSave+1).zfill(3) + 'of' + str(cfg.recon.sliceCount).zfill(3)
        fileName = cfg.resultsName + '_' + sliceNumberString + '_' + imageVolume3D_size_string + '.raw'
        sliceToSave = imageVolume3D[:, :, sliceIndexToSave]
        sliceToSave = sliceToSave.copy(order='C')
        rawwrite(fileName, sliceToSave)


def displayImagePictures(cfg, imageVolume3D):

    cfg = drawImages('screen', cfg, imageVolume3D)
    
    return cfg


def saveImagePictureFiles(cfg, imageVolume3D):

    print('* Saving the recon results to individual .png files...')

    cfg = drawImages('file', cfg, imageVolume3D)
    
    return cfg

def drawImages(drawTo, cfg, imageVolume3D):

    sliceIndicesToDraw = range(0, cfg.recon.sliceCount)

    if hasattr(cfg, 'displayWindowMin') and hasattr(cfg, 'displayWindowMax'):
        # If displayWindowMin and displayWindowMax are passed in, use them.
        displayWindowMin = cfg.displayWindowMin
        displayWindowMax = cfg.displayWindowMax
    else:
        # Otherwise, find displayWindowMin and displayWindowMax.
        displayWindowMin = np.min(imageVolume3D)
        displayWindowMax = np.max(imageVolume3D)

    displayWindow = displayWindowMax - displayWindowMin
    displayLevel = displayWindow/2
      
    for sliceIndexToDraw in sliceIndicesToDraw:
        sliceToDraw = imageVolume3D[:, :, sliceIndexToDraw]
        sliceToDraw = sliceToDraw.copy(order='C')
        sliceNumberString = 'slice' + str(sliceIndexToDraw+1).zfill(3) + 'of' + str(cfg.recon.sliceCount).zfill(3)
        fileName = cfg.resultsName + '_' + sliceNumberString + '.png'
        plt.figure(int(sliceIndexToDraw+1))
        plt.imshow(sliceToDraw, cmap='gray', vmin=displayWindowMin, vmax=displayWindowMax)
        sliceString = "slice " + str(sliceIndexToDraw+1) + " of " + str(cfg.recon.sliceCount) + "\n"
        if cfg.recon.unit == 'HU':
            formatString = "W/L = {}/{} {:s}; cfg.physics.monochromatic = {};"
        if cfg.recon.unit == '/cm':
            formatString = "W/L = {}/{} {:s}; cfg.physics.monochromatic = {};"
        if cfg.recon.unit == '/mm':
            formatString = "W/L = {}/{} {:s}; cfg.physics.monochromatic = {};"
        string1 = formatString.format(displayWindow, displayLevel, cfg.recon.unit, cfg.physics.monochromatic)
        string2 = "cfg.physics.enableElectronicNoise = {}; cfg.protocol.spectrumScaling = {};".format(cfg.physics.enableElectronicNoise, cfg.protocol.spectrumScaling)
        string3 = "cfg.physics.enableQuantumNoise = {}; cfg.protocol.mA = {}".format(cfg.physics.enableQuantumNoise, cfg.protocol.mA)
        plt.title(string1 + "\n" + string2 + "\n" + string3, fontsize=10)

        if drawTo == 'file':
            plt.savefig(fileName, bbox_inches='tight')
            plt.close()
        elif drawTo == 'screen':
            plt.draw()
            
    if drawTo == 'screen':
        plt.pause(1)
        if cfg.waitForKeypress:
            print('********************************************')
            print('* Press Enter to close images and continue *')
            input('********************************************')
        plt.close('all')

    cfg.displayWindowMin = displayWindowMin
    cfg.displayWindowMax = displayWindowMax
    
    return cfg
