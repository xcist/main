from gecatsim.pyfiles.CommonTools import *
import unittest.mock
from gecatsim.pyfiles import CommonTools
from gecatsim.pyfiles.SetFocalspot import *

class Test_SetFocalSpot(unittest.TestCase):

    def test_SetFocalSpot(self):
        cfg = CommonTools.CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                              "../examples/cfg/Protocol_Sample_axial")

        cfg = SetFocalspot(cfg)
        assert(cfg is not None)

    def test_GetDefaultWidthLength(self):

        shape = "uniform"
        width, length = GetDefaultWidthLength(shape)
        assert(width==1)
        assert(length==1)

        shape = "Gaussian"
        width, length = GetDefaultWidthLength(shape)
        assert(width==1)
        assert(length==1)

    def test_GetIntensity(self):

        cfg = CommonTools.CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                              "../examples/cfg/Protocol_Sample_axial")
        
        weights = GetIntensity(cfg)
        assert(weights is not None)
        assert(len(weights)>0)
