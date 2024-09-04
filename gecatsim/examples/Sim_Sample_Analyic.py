# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

###------------ import XCIST-CatSim
import gecatsim as xc

##--------- Initialize 
ct = xc.CatSim("./cfg/Phantom_Sample_Analytic")  # initialization

##--------- Make changes to parameters (optional)
# ct.phantom.filename = 'water20.ppm'
ct.phantom.filename = 'CTDI_16cm_WaterAirPEBoneChambers.ppm'

ct.resultsName = "test_Analytic"
ct.protocol.viewsPerRotation = 500
ct.protocol.viewCount = ct.protocol.viewsPerRotation
ct.protocol.stopViewId = ct.protocol.viewCount-1

ct.physics.enableQuantumNoise = 0
ct.physics.enableElectronicNoise = 0

ct.physics.callback_post_log = 'Prep_BHC_Accurate'
ct.physics.EffectiveMu = 0.2
ct.physics.BHC_poly_order = 5
ct.physics.BHC_max_length_mm = 220
ct.physics.BHC_length_step_mm = 10

ct.physics.colSampleCount = 2
ct.physics.rowSampleCount = 2
ct.physics.srcXSampleCount = 2
ct.physics.srcYSampleCount = 2
ct.physics.viewSampleCount = 1

##--------- Run simulation
ct.run_all()  # run the scans defined by protocol.scanTypes


##--------- Show results
import numpy as np
import matplotlib.pyplot as plt

prep = xc.rawread(ct.resultsName+'.prep', [ct.protocol.viewCount, ct.scanner.detectorRowCount, ct.scanner.detectorColCount], 'float')
prep = prep[-1, :, :]
plt.plot(prep[1, :])
plt.show()
