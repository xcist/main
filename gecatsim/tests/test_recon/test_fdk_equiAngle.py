import gecatsim as xc
from gecatsim.pyfiles.CommonTools import *
from gecatsim.reconstruction.pyfiles.recon import *
from gecatsim.reconstruction.pyfiles.fdk_equiAngle import *
import unittest.mock

class Test_fdk_equiAngle(unittest.TestCase):
    def test_fdk_equiAngle_1(self):

        ##--------- Initialize
        cfg = xc.CatSim("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic", "../examples/cfg/Protocol_Sample_axial")

        ##--------- Make changes to parameters (optional)
        cfg.resultsName = "test"
        cfg.protocol.viewsPerRotation = 500 
        cfg.protocol.viewCount = cfg.protocol.viewsPerRotation
        cfg.protocol.stopViewId = cfg.protocol.viewCount-1
        # cfg.protocol.scanTypes = [1, 0, 0]  # flags for airscan, offset scan, phantom scan
        # cfg.load_cfg("Protocol_Sample_axial", "Physics_Sample", "Recon_Sample_2d")  # new cfg overrides existing parameters

        cfg.protocol.mA = 800 
        cfg.scanner.detectorRowsPerMod = 4 
        cfg.scanner.detectorRowCount = cfg.scanner.detectorRowsPerMod

        cfg.recon.fov = 300.0
        cfg.recon.sliceCount = 4        # number of slices to reconstruct
        cfg.recon.sliceThickness = 0.568  # reconstruction inter-slice interval (in mm)

        ##--------- Run simulation
        cfg.run_all()  # run the scans defined by protocol.scanTypes

        ##--------- Reconstruction
        cfg.do_Recon = 1 
        prep = load_prep(cfg)

        # cfg.recon.reconType is the recon function's name
        imageVolume3D = feval("gecatsim.reconstruction.pyfiles." + cfg.recon.reconType, cfg, prep)
        assert(sys.getsizeof(imageVolume3D) == 144)

        imageVolume3D = scaleReconData(cfg, imageVolume3D)
        assert(sys.getsizeof(imageVolume3D) == 4194448)

