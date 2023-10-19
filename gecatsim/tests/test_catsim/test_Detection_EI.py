from gecatsim.pyfiles.Detection_EI import Detection_EI
import unittest.mock
from gecatsim.pyfiles.CommonTools import *
import numpy as np

class test_Detection_EI(unittest.TestCase):

    def test_Detection_EI(self, signal_mock):
        cfg = CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                    "../examples/cfg/Protocol_Sample_axial")

        cfg.sim.startViewId = 0
        cfg.sim.stopViewId = 2
        cfg.sim.enableQuantumNoise = 1
        cfg.sim.subViewCount = 1
        cfg.det.totalNumCells = 5
        cfg.det.cosBetas = np.ones([cfg.det.totalNumCells, 1], dtype=np.single)

        cfg = feval(cfg.protocol.spectrumCallback, cfg)

        cfg.thisSubView = np.float32(np.random.random([cfg.det.totalNumCells, cfg.sim.Evec.size]) * 100)
        cfg.thisSubView[1:4, 3:5] = 0
        viewId = 0
        subViewId = 0
        cfg.sim.eNoise = 3500

        assert cfg.sim.Evec is not None
        assert cfg.thisSubView is not None

        cfg = Detection_EI(cfg, viewId, subViewId)
        assert cfg.sim.Wvec is not None
        assert cfg.thisView is not None
