import catsim as xc
from reconstruction.pyfiles import recon

ct = xc.CatSim("all_in_one")  # initialization

ct.resultsName = "out"

ct.run_all()  # run the scans defined by protocol.scanTypes
# remove pretrigger 100 views
oldprep = xc.rawread('out.prep', [ct.protocol.viewCount, ct.scanner.detectorRowCount, ct.scanner.detectorColCount], 'float')
xc.rawwrite('out.prep', oldprep[100:])

if ct.physics.monochromatic>0:
        ct.recon.mu = xc.GetMu('water', ct.physics.monochromatic)[0]/10

cfg = ct.get_current_cfg();
cfg.recon.do_Recon = 1
cfg.waitForKeypress = 0
#recon.recon(cfg)
