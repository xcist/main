# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

###------------ import XCIST-CatSim
import numpy as np
import gecatsim as xc
import gecatsim.reconstruction.pyfiles.recon as recon

##--------- Initialize
ct = xc.CatSim("./cfg/Phantom_Sample_Polygonal")  # initialization

##--------- Make changes to parameters (optional)
ct.resultsName = "test_Polygonal"

ct.phantom.filename = '../phantom/female_adult_average_lung_lesions_reduced.nrb'

# Here the phantom position and scaling are only for illustrating purpose due to the object's small size.
ct.phantom.centerOffset = [54, 51, -102]      # offset of phantom center relative to origin (in mm)
ct.phantom.scale = 3
ct.phantom.centerOffset = np.array(ct.phantom.centerOffset)*ct.phantom.scale

ct.protocol.viewsPerRotation = 500
ct.protocol.viewCount = ct.protocol.viewsPerRotation
ct.protocol.stopViewId = ct.protocol.viewCount-1

ct.physics.colSampleCount = 2
ct.physics.rowSampleCount = 1
ct.physics.srcXSampleCount = 2
ct.physics.srcYSampleCount = 1
ct.physics.viewSampleCount = 1

##--------- Run simulation
ct.run_all()  # run the scans defined by protocol.scanTypes

##--------- Reconstruction
ct.recon.fov = 50.0
ct.recon.sliceCount = ct.scanner.detectorRowCount   # number of slices to reconstruct
ct.recon.sliceThickness = 0.568  # reconstruction inter-slice interval (in mm)
ct.do_Recon = 1
recon.recon(ct)

##--------- Show results
import matplotlib.pyplot as plt

imgFname = "%s_%dx%dx%d.raw" %(ct.resultsName, ct.recon.imageSize, ct.recon.imageSize, ct.recon.sliceCount)
img = xc.rawread(imgFname, [ct.recon.sliceCount, ct.recon.imageSize, ct.recon.imageSize], 'float')
plt.imshow(img[np.int32(ct.recon.sliceCount/2),:,:], cmap='gray', vmin=-1000, vmax=200)
plt.show()
