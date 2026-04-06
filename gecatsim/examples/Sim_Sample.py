# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

###------------ import XCIST-CatSim
# import the package from the current directory (instead of the installed package) to ensure the latest code is used
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))  # add the parent directory of the current file to the search path
import gecatsim as xc


##--------- Initialize
#my_path = xc.pyfiles.CommonTools.my_path
# add any additional search directories
#my_path.add_search_path("my-experiments")

example_dir = os.path.dirname(os.path.abspath(__file__))
ct = xc.CatSim(
	os.path.join(example_dir, "cfg", "Phantom_Sample"),
	os.path.join(example_dir, "cfg", "Scanner_Sample_generic"),
	os.path.join(example_dir, "cfg", "Protocol_Sample_axial"),
)  # initialization

##--------- Make changes to parameters (optional)
ct.resultsName = "test"
ct.protocol.viewsPerRotation = 50
ct.protocol.viewCount = ct.protocol.viewsPerRotation
ct.protocol.stopViewId = ct.protocol.viewCount-1
# ct.protocol.scanTypes = [1, 0, 0, 0]  # flags for airscan, offset scan, phantom scan, prep
# ct.load_cfg("Protocol_Sample_axial", "Physics_Sample", "Recon_Sample_2d")  # new cfg overrides existing parameters

##--------- Run simulation
ct.run_all()  # run the scans defined by protocol.scanTypes


##--------- Show results
import matplotlib.pyplot as plt

prep = xc.rawread(ct.resultsName+'.prep', [ct.protocol.viewCount, ct.scanner.detectorRowCount, ct.scanner.detectorColCount], 'float')
prep = prep[-1, :, :]
plt.plot(prep[7, :])
plt.show()
