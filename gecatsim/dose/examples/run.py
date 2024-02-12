import sys, os
import gecatsim as xc
import gecatsim.reconstruction.pyfiles.recon as recon
from gecatsim.dose.pyfiles.catdoserecon import catdoserecon


doSim = True
doRecon = True
doDose = True

ct = xc.CatSim("dose_physics", "dose_phantom", "dose_recon", "dose_scanner", "dose_protocol", "dose_doserecon")  # initialization

# add any additional search directories
my_path = xc.pyfiles.CommonTools.my_path
my_path.add_search_path("my_data")

ct.resultsName = "out" # for sim and recon

if doSim:
    ct.run_all()  # run the scans defined by protocol.scanTypes

if doRecon:
    ct.do_Recon = 1
    recon.recon(ct)

if doDose:
    ct.dose.imageFileName = 'out_512x512x4.raw'    # name of recon image
    ct.dose.doseFileName = 'dose'                  # for doserecon
    cfg = ct.get_current_cfg()
    dosevol = catdoserecon(cfg=cfg)

## check results
import matplotlib.pyplot as plt
doseVolFname = "%s.mGy" %(ct.dose.doseFileName)
img = xc.rawread(doseVolFname, [ct.recon.sliceCount, ct.dose.nVoxel, ct.dose.nVoxel], 'float')
plt.imshow(img[0,:,:])
plt.show()
