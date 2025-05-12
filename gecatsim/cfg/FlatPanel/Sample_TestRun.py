# %%
# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

###------------ import XCIST-CatSim
import gecatsim as xc
import numpy as np
import matplotlib.pyplot as plt

##--------- Initialize 
ct = xc.CatSim("Phantom_Sample_FlatPanel", "Scanner_Sample_FlatPanel", "Physics_Sample_FlatPanel", "Protocol_Sample_FlatPanel")  # initialization
# %%

##--------- Make changes to parameters (optional)
# Name of the outputs
ct.resultsName = "airscan_test"

# ct.protocol.viewsPerRotation = 50
# ct.protocol.viewCount = ct.protocol.viewsPerRotation
# ct.protocol.stopViewId = ct.protocol.viewCount-1
# ct.protocol.scanTypes = [1, 0, 0]  # flags for airscan, offset scan, phantom scan
#ct.load_cfg("Protocol_Sample_axial", "Physics_Sample", "Phantom_Sample", "Scanner_Sample_generic", "Recon_Sample_2d")  # new cfg overrides existing parameters

##--------- Run simulation
ct.run_all()  # run the scans defined by protocol.scanTypes
# %%

##--------- Show results
airscan = xc.rawread(ct.resultsName+'.air', [500, 500], 'float')
a = airscan[250,:]
plt.figure()
plt.plot(a)
plt.show()

# %%
