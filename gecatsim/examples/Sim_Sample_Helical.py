# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

###------------ import XCIST-CatSim
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
import gecatsim as xc
from gecatsim.reconstruction.pyfiles import recon

##--------- Initialize 
example_dir = os.path.dirname(os.path.abspath(__file__))
ct = xc.CatSim(os.path.join(example_dir, "cfg", "Phantom_Sample_Analytic"),
                 os.path.join(example_dir, "cfg", "Protocol_Sample_Helical"),
                 os.path.join(example_dir, "cfg", "Scanner_Sample_generic"),
                 os.path.join(example_dir, "cfg", "Physics_Sample"),
                 os.path.join(example_dir, "cfg", "Recon_Sample_Helical"),

        )  # initialization

##--------- Make changes to parameters (optional)
# ct.phantom.filename = 'water20.ppm'
ct.phantom.filename = 'CTDI_16cm_WaterAirPEBoneChambers.ppm'

ct.resultsName = "test_Helical"

##--------- Run simulation
ct.run_all()  # run the scans defined by protocol.scanTypes


if ct.physics.monochromatic>0:
        ct.recon.mu = xc.GetMu('water', ct.physics.monochromatic)[0]/10

cfg = ct.get_current_cfg();
cfg.do_Recon = 1
cfg.waitForKeypress = 0
recon.recon(cfg)
