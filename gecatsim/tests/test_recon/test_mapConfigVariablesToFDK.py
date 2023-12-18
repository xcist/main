import gecatsim as xc
from gecatsim.pyfiles.CommonTools import *
from gecatsim.reconstruction.pyfiles.recon import *
from gecatsim.reconstruction.pyfiles.mapConfigVariablesToFDK import *
import unittest.mock

class Test_mapConfigVariablesToFDK(unittest.TestCase):
    def test_mapConfigVariablesToFDK_1(self):
        length = 4

        cfg = xc.CatSim("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic", "../examples/cfg/Protocol_Sample_axial")

        sid, sdd, nMod, rowSize, modWidth, dectorYoffset, dectorZoffset, \
        fov, imageSize, sliceCount, sliceThickness, centerOffset, startView, rotdir, kernelType \
            = mapConfigVariablesToFDK(cfg)

        assert(fov==500)
        assert(imageSize==512)

        assert(1==1)
