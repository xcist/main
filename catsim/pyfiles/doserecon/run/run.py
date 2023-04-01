import catsim as xc
from reconstruction.pyfiles import recon
import sys
import os

doSim = True
doRecon = True
doDose = True

ct = xc.CatSim("dose_physics", "dose_phantom", "dose_recon", "dose_scanner", "dose_protocol", "dose_doserecon")  # initialization
ct.resultsName = "out" # for sim and recon

if doSim:
    ct.run_all()  # run the scans defined by protocol.scanTypes

if doRecon:
    cfg = ct.get_current_cfg();
    cfg.do_Recon = 1
    cfg.waitForKeypress = 0
    recon.recon(cfg)

if doDose:
    cfg = ct.get_current_cfg();
    cfg.sim.isOffsetScan = 0
    cfg.results_basename = "doserecon" # for doserecon
    cfg.imageFileName = 'out';  #% name of recon image
    #ct.doseFileName = 'doserecon';
    from catsim.pyfiles.doserecon.pyfiles.catdoserecon import catdoserecon
    dosevol = catdoserecon(cfg=cfg);
