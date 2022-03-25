# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

from catsim.CommonTools import *
import matplotlib.pyplot as plt
# Need to import new recons as they are added
from reconstruction.fdk_equiAngle import fdk_equiAngle


def recon(cfg):

    prep = load_prep(cfg)

    # The following line doesn't work - need to fix it when new recons are added.
    # imageVolume3D = feval('reconstruction.' + cfg.recon.reconType, cfg, prep)

    # A hack until the previos line is fixed.
    imageVolume3D = fdk_equiAngle(cfg, prep)

    if cfg.recon.saveImageVolume or cfg.recon.saveSingleImages:
        imageVolume3D = scaleReconData(cfg, imageVolume3D)
        if cfg.recon.saveImageVolume:
            saveImageVolume(cfg, imageVolume3D)
        if cfg.recon.saveSingleImages:
            saveSingleImages(cfg, imageVolume3D)

    if cfg.recon.displayImagePictures:
        displayImagePictures(cfg, imageVolume3D)

    if cfg.recon.saveImagePictureFiles:
        saveImagePictureFiles(cfg, imageVolume3D)

     
def load_prep(cfg):

    print('* Loading the projection data...')
    prep = rawread(cfg.resultsName + '.prep',
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
        raise Exception('******** Error! An recon unit was specified: {:s}. ********'.format(cfg.recon.unit))


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

    sliceIndicesToDisplay = range(0, cfg.recon.sliceCount)
    for sliceIndexToDisplay in sliceIndicesToDisplay:
        sliceToDisplay = imageVolume3D[:, :, sliceIndexToDisplay]
        sliceToDisplay = sliceToDisplay.copy(order='C')
        plt.figure(int(sliceIndexToDisplay+1))
        plt.imshow(sliceToDisplay, cmap='gray')
        plt.title("slice " + str(sliceIndexToDisplay+1) + " of " + str(cfg.recon.sliceCount))

    plt.draw()
    plt.pause(1)
    print('********************************************')
    print('* Press Enter to close images and continue *')
    input('********************************************')
    plt.close('all')


def saveImagePictureFiles(cfg, imageVolume3D):

    print('* Saving the recon results to individual .png files...')

    sliceIndicesToSave = range(0, cfg.recon.sliceCount)
    for sliceIndexToSave in sliceIndicesToSave:
        sliceToSave = imageVolume3D[:, :, sliceIndexToSave]
        sliceToSave = sliceToSave.copy(order='C')
        sliceNumberString = 'slice' + str(sliceIndexToSave+1).zfill(3) + 'of' + str(cfg.recon.sliceCount).zfill(3)
        fileName = cfg.resultsName + '_' + sliceNumberString + '.png'
        plt.figure  
        plt.imshow(sliceToSave, cmap='gray')
        plt.title("slice " + str(sliceIndexToSave+1) + " of " + str(cfg.recon.sliceCount))
        plt.savefig(fileName, bbox_inches='tight')
        plt.close()
