# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

import numpy as np
import matplotlib.pyplot as plt

###------------ import XCIST-CatSim
import catsim as xc
import recon


##--------- Initialize 
ct = xc.CatSim('../cfg/Physics_VCT','../cfg/Protocol_VCT','../cfg/Scanner_VCT')  # initialization

##--------- Make changes to parameters (optional)
ct.resultsName = "../sim/View/test"

ct.physics.enableQuantumNoise = 1
ct.physics.enableElectronicNoise = 1
ct.protocol.mA = 800

ct.protocol.spectrumFilename = "tungsten_tar7_120_unfilt.dat"
ct.physics.energyCount = 12
# ct.physics.monochromatic = 70

ct.phantom.filename = 'C:/Users/302017896/Downloads/Duke_phantom/CatSim_Adult_Male_Standard_Lung_Phantom_1700x1050x900/Slab_400/reduced/adult_male_standard_lung.json'

##--------- Run simulation
# cfg = ct.run_all()  # run the scans defined by protocol.scanTypes
cfg = ct.get_current_cfg();
recon.FDK(cfg)


##--------- Show results
# import numpy as np
# import matplotlib.pyplot as plt

# prep = xc.rawread(ct.resultsName+'.prep', [ct.protocol.viewCount, ct.scanner.detectorColCount, ct.scanner.detectorRowCount], 'float')
# prep = prep[-1, :, :]
# plt.plot(prep[:, 7])
# plt.show()
