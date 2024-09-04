# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

###------------ import XCIST-CatSim
import gecatsim as xc


##--------- Initialize 
ct = xc.CatSim()  # initialization

##--------- Make changes to parameters (optional)
ct.resultsName = "test_sc"
ct.physics.scatterCallback = "Scatter_ConvolutionModel"  # scatter model
ct.physics.scatterKernelCallback = ""  # using default
ct.physics.scatterScaleFactor = 1

ct.physics.callback_pre_log = "Scatter_Correction"  # scatter correction
ct.physics.scatterCorrectionScaleFactor = 1
ct.physics.scatterCorrectionSaveView = 1  # save the corrected phantom scan

ct.protocol.viewsPerRotation = 1
ct.protocol.viewCount = ct.protocol.viewsPerRotation
ct.protocol.stopViewId = ct.protocol.viewCount-1

ct.scanner.detectorRowsPerMod = 64
ct.scanner.detectorRowCount = ct.scanner.detectorRowsPerMod

##--------- Run simulation
ct.run_all()  # run the scans defined by protocol.scanTypes


##--------- Show results
import numpy as np
import matplotlib.pyplot as plt

prep = xc.rawread(ct.resultsName+'.prep', [ct.protocol.viewCount, ct.scanner.detectorRowCount, ct.scanner.detectorColCount], 'float')
prep = prep[-1, :, :]
plt.plot(prep[7, :])
plt.show()
