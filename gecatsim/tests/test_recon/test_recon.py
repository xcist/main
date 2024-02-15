from gecatsim.pyfiles.CommonTools import *
import gecatsim.reconstruction.pyfiles.recon as recon
import gecatsim as xc
import unittest.mock
import sys

class Test_recon(unittest.TestCase):
    def test_recon_1(self):

        ##--------- Initialize
        ct = xc.CatSim("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic", "../examples/cfg/Protocol_Sample_axial")

        ##--------- Make changes to parameters (optional)
        ct.resultsName = "test"
        ct.protocol.viewsPerRotation = 500
        ct.protocol.viewCount = ct.protocol.viewsPerRotation
        ct.protocol.stopViewId = ct.protocol.viewCount-1
        # ct.protocol.scanTypes = [1, 0, 0]  # flags for airscan, offset scan, phantom scan
        # ct.load_cfg("Protocol_Sample_axial", "Physics_Sample", "Recon_Sample_2d")  # new cfg overrides existing parameters

        ct.protocol.mA = 800
        ct.scanner.detectorRowsPerMod = 4
        ct.scanner.detectorRowCount = ct.scanner.detectorRowsPerMod

        ct.recon.fov = 300.0
        ct.recon.sliceCount = 4        # number of slices to reconstruct
        ct.recon.sliceThickness = 0.568  # reconstruction inter-slice interval (in mm)

        ##--------- Run simulation
        #ct.run_all()  # run the scans defined by protocol.scanTypes

        ##--------- Reconstruction
        ct.do_Recon = 1
        retval = recon.recon(ct)

        assert(sys.getsizeof(retval)>0)
        assert(ct==retval)