# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

###------------ import XCIST-CatSim
import gecatsim as xc
import gecatsim.reconstruction.pyfiles.recon as recon
from gecatsim.pyfiles.CommonTools import *

##--------- Initialize
# add any additional search directories
#my_path = xc.pyfiles.CommonTools.my_path
#my_path.add_search_path("my-experiments")

ct = xc.CatSim("./cfg/Scanner_PCCT", "./cfg/Phantom_Sample", "./cfg/Protocol_Sample_axial",'./cfg/Physics_Sample')  # initialization

##--------- PCCT
ct.scanner.detectorMaterial = "CZT"            # detector sensor material
ct.scanner.detectorDepth = 1.6                 # detector sensor depth (in mm)

ct.scanner.detectionCallback = "Detection_PC"  # name of function that defines the detection process (conversion from X-rays to detector signal)
ct.scanner.detectionResponseFilename = 'PC_spectral_response_CZT0.25x0.25x1.6.mat'   # name of the response data file
ct.scanner.detectorBinThreshold = [20, 30, 40, 60, 80, 100, 160]  # energy thresholds (keV), n bins has n+1 thresholds; the first and last are the min and max energy thresholds.
ct.scanner.detectorSumBins = 0                 # 1: sum all bins (gray scale output), data dim [view row col]; 0: output multiple bins [view row col bin]

##--------- Make changes to parameters (optional)
ct.resultsName = "test"
ct.protocol.viewsPerRotation = 500
ct.protocol.viewCount = ct.protocol.viewsPerRotation
ct.protocol.stopViewId = ct.protocol.viewCount-1

ct.protocol.mA = 500
ct.protocol.spectrumFilename = "tungsten_tar7.0_120_filt.dat"
ct.physics.energyCount = 120

ct.protocol.bowtie = "large.txt"
ct.protocol.flatFilter = ['Al', 3.0]

ct.scanner.detectorRowsPerMod = 2
ct.scanner.detectorRowCount = ct.scanner.detectorRowsPerMod

ct.recon.fov = 300.0
ct.recon.sliceCount = 1        # number of slices to reconstruct
ct.recon.sliceThickness = 0.1421  # reconstruction inter-slice interval (in mm)

##--------- Run simulation
if not ct.scanner.detectorSumBins:
    ct.do_prep = 0
ct.run_all()  # run the scans defined by protocol.scanTypes

##--------- Reconstruction
if ct.scanner.detectorSumBins==1:
    ct.do_Recon = 1
    recon.recon(ct)


##--------- Show results
import matplotlib.pyplot as plt

if ct.scanner.detectorSumBins==1:
    imgFname = "%s_%dx%dx%d.raw" %(ct.resultsName, ct.recon.imageSize, ct.recon.imageSize, ct.recon.sliceCount)
    img = xc.rawread(imgFname, [ct.recon.sliceCount, ct.recon.imageSize, ct.recon.imageSize], 'float')
    plt.imshow(img[0,:,:], cmap='gray', vmin=-200, vmax=200)
    plt.show()
else:
    scanFname = "%s.air" % ct.resultsName
    nBin = len(ct.scanner.detectorBinThreshold)-1
    air = xc.rawread(scanFname, [ct.scanner.detectorRowCount, ct.scanner.detectorColCount, nBin], 'float')
    plt.plot(air[0, :, :])
    plt.legend(np.arange(1,nBin+1))
    plt.show()

    plt.plot(air[0,int(ct.scanner.detectorColCount/2)-1,:])
    plt.show()
