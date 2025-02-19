# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import numpy as np
import matplotlib.pyplot as plt

###------------ import XCIST-CatSim
import gecatsim as xc
import gecatsim.reconstruction.pyfiles.recon as recon

##--------- Initialize
#my_path = xc.pyfiles.CommonTools.my_path
# add any additional search directories
#my_path.add_search_path("my-experiments")

ct = xc.CatSim("./cfg/Phantom_Sample", "./cfg/Scanner_PCCT", "./cfg/Physics_Sample")  # initialization

##--------- Make changes to parameters (optional)
ct.resultsName = "test"
ct.protocol.viewsPerRotation = 500
ct.protocol.viewCount = ct.protocol.viewsPerRotation
ct.protocol.stopViewId = ct.protocol.viewCount-1
ct.protocol.mA = 800

ct.scanner.detectorColSize = 1                 # detector column pitch or size (in mm)
ct.scanner.detectorRowSize = 1                 # detector row pitch or size (in mm)
ct.scanner.detectorColCount = 900              # total number of detector columns
ct.scanner.detectorRowsPerMod = 4
ct.scanner.detectorRowCount = ct.scanner.detectorRowsPerMod

ct.recon.fov = 300.0
ct.recon.sliceCount = 4        # number of slices to reconstruct
ct.recon.sliceThickness = 0.568  # reconstruction inter-slice interval (in mm)

##--------- Run simulation
if not ct.scanner.detectorSumBins:
    ct.do_prep = 0
ct.run_all()  # run the scans defined by protocol.scanTypes

##--------- Prep and recon for each bin
nBin = len(ct.scanner.detectorBinThreshold)-1
nCol = ct.scanner.detectorColCount
nRow = ct.scanner.detectorRowCount
nView = ct.protocol.viewCount

airscan = xc.rawread("%s.air" % ct.resultsName, [nRow, nCol, nBin], 'float')
offsetscan = xc.rawread("%s.offset" % ct.resultsName, [nRow, nCol, nBin], 'float')
phantomscan = xc.rawread("%s.scan" % ct.resultsName, [nView, nRow, nCol, nBin], 'float')
airscan -= offsetscan
phantomscan -= offsetscan

for binId in range(nBin):
    prep = np.where(phantomscan[:,:,:,binId]==0,0,-np.log(phantomscan[:,:,:,binId]/airscan[:,:,binId]))
    prep = np.where(prep<0,0,prep)
    prep = np.where(prep>20,20,prep)

    newFname = "test_bin%s" % binId
    xc.rawwrite(newFname+'.prep', prep)
    ct.resultsName = newFname

    ##--------- Reconstruction
    ct.do_Recon = 1
    recon.recon(ct)


##--------- Show results
imgFname = "%s_%dx%dx%d.raw" %(ct.resultsName, ct.recon.imageSize, ct.recon.imageSize, ct.recon.sliceCount)
img = xc.rawread(imgFname, [ct.recon.sliceCount, ct.recon.imageSize, ct.recon.imageSize], 'float')
plt.imshow(img[2,:,:], cmap='gray', vmin=-200, vmax=200)
plt.show()
