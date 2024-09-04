# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

###------------ import XCIST-CatSim
import gecatsim as xc
import gecatsim.reconstruction.pyfiles.recon as recon


##--------- Initialize
#my_path = xc.pyfiles.CommonTools.my_path
# add any additional search directories
#my_path.add_search_path("my-experiments")

ct = xc.CatSim("./cfg/Phantom_Sample", "./cfg/Scanner_Sample_generic", "./cfg/Protocol_Sample_axial")  # initialization

##--------- Make changes to parameters (optional)
ct.resultsName = "test"
ct.protocol.viewsPerRotation = 500
ct.protocol.viewCount = ct.protocol.viewsPerRotation
ct.protocol.stopViewId = ct.protocol.viewCount-1
# ct.protocol.scanTypes = [1, 0, 0]  # flags for airscan, offset scan, phantom scan
# ct.load_cfg("Protocol_Sample_axial", "Physics_Sample", "Recon_Sample_2d")  # new cfg overrides existing parameters

ct.protocol.mA = 800
ct.scanner.detectorRowsPerMod = 4
ct.scanner.detectorRowCount = ct.scanner.detectorRowsPerMod

ct.recon.fov = 300.0
ct.recon.sliceCount = 4        # number of slices to reconstruct
ct.recon.sliceThickness = 0.568  # reconstruction inter-slice interval (in mm)

##--------- Run simulation
ct.run_all()  # run the scans defined by protocol.scanTypes

##--------- Reconstruction
ct.do_Recon = 1
recon.recon(ct)


##--------- Show results
import matplotlib.pyplot as plt

imgFname = "%s_%dx%dx%d.raw" %(ct.resultsName, ct.recon.imageSize, ct.recon.imageSize, ct.recon.sliceCount)
img = xc.rawread(imgFname, [ct.recon.sliceCount, ct.recon.imageSize, ct.recon.imageSize], 'float')
plt.imshow(img[2,:,:], cmap='gray', vmin=-200, vmax=200)
plt.show()
