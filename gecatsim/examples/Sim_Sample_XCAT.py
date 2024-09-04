# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

###------------ import XCIST-CatSim
import numpy as np
import gecatsim as xc
import gecatsim.reconstruction.pyfiles.recon as recon

##--------- Initialize 
ct = xc.CatSim("./cfg/Phantom_Sample_XCAT")  # initialization

##--------- Make changes to parameters (optional)
ct.resultsName = "test_XCAT"

ct.phantom.filename = '../phantom/vmale50_chest_less_surfaces.nrb'

ct.protocol.viewsPerRotation = 500
ct.protocol.viewCount = ct.protocol.viewsPerRotation
ct.protocol.stopViewId = ct.protocol.viewCount-1

# ct.physics.enableQuantumNoise = 0
# ct.physics.enableElectronicNoise = 0

ct.physics.colSampleCount = 2
ct.physics.rowSampleCount = 1
ct.physics.srcXSampleCount = 2
ct.physics.srcYSampleCount = 1
ct.physics.viewSampleCount = 1

##--------- Run simulation
ct.run_all()  # run the scans defined by protocol.scanTypes

##--------- Reconstruction
ct.recon.fov = 400.0
ct.recon.sliceCount = ct.scanner.detectorRowCount   # number of slices to reconstruct
ct.recon.sliceThickness = 0.568  # reconstruction inter-slice interval (in mm)
ct.do_Recon = 1
recon.recon(ct)

##--------- Show results
import matplotlib.pyplot as plt

# prep = xc.rawread(ct.resultsName+'.prep', [ct.protocol.viewCount, ct.scanner.detectorRowCount, ct.scanner.detectorColCount], 'float')
# prep = prep[-1, :, :]
# plt.plot(prep[1, :])
# plt.show()

imgFname = "%s_%dx%dx%d.raw" %(ct.resultsName, ct.recon.imageSize, ct.recon.imageSize, ct.recon.sliceCount)
img = xc.rawread(imgFname, [ct.recon.sliceCount, ct.recon.imageSize, ct.recon.imageSize], 'float')
plt.imshow(img[np.int32(ct.recon.sliceCount/2),:,:], cmap='gray', vmin=-1000, vmax=200)
plt.show()
