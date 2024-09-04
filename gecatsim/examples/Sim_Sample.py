# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

###------------ import XCIST-CatSim
import gecatsim as xc


##--------- Initialize
#my_path = xc.pyfiles.CommonTools.my_path
# add any additional search directories
#my_path.add_search_path("my-experiments")

ct = xc.CatSim("./cfg/Phantom_Sample", "./cfg/Scanner_Sample_generic", "./cfg/Protocol_Sample_axial")  # initialization

##--------- Make changes to parameters (optional)
ct.resultsName = "test"
ct.protocol.viewsPerRotation = 50
ct.protocol.viewCount = ct.protocol.viewsPerRotation
ct.protocol.stopViewId = ct.protocol.viewCount-1
# ct.protocol.scanTypes = [1, 0, 0]  # flags for airscan, offset scan, phantom scan
# ct.load_cfg("Protocol_Sample_axial", "Physics_Sample", "Recon_Sample_2d")  # new cfg overrides existing parameters

##--------- Run simulation
ct.run_all()  # run the scans defined by protocol.scanTypes


##--------- Show results
import matplotlib.pyplot as plt

prep = xc.rawread(ct.resultsName+'.prep', [ct.protocol.viewCount, ct.scanner.detectorRowCount, ct.scanner.detectorColCount], 'float')
prep = prep[-1, :, :]
plt.plot(prep[7, :])
plt.show()
