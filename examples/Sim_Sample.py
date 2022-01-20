# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

###------------ import XCIST-CatSim
import catsim as xc


##--------- Initialize 
ct = xc.CatSim("Phantom_Sample", "Scanner_Sample_generic", "Protocol_Sample_axial")  # initialization

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
import numpy as np
import matplotlib.pyplot as plt

prep = xc.rawread(ct.resultsName+'.prep', [ct.protocol.viewCount, ct.scanner.detectorRowCount, ct.scanner.detectorColCount], 'float')
prep = prep[-1, :, :]
plt.plot(prep[7, :])
plt.show()
