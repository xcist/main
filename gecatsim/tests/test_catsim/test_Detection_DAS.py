import unittest.mock
from gecatsim.pyfiles.CommonTools import *
from gecatsim.pyfiles.Detection_DAS import Detection_DAS


class TestDetection_DAS(unittest.TestCase):

    def test_detection_DAS(self):
        cfg = CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                  "../examples/cfg/Protocol_Sample_axial")

        cfg.sim.eNoise = 3500.0

        viewIn = np.array([1])
        viewId = 1

        viewOut = Detection_DAS(viewIn, viewId, cfg)

        assert viewOut is not None

