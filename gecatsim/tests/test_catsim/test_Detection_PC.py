import unittest.mock
from unittest.mock import patch

from gecatsim.pyfiles.CommonTools import *
from gecatsim.pyfiles.Detection_PC import Detection_PC


class TestDetection_PC(unittest.TestCase):

    def test_Detection_PC(self):
        cfg = CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                  "../examples/cfg/Protocol_Sample_axial")

        cfg.scanner.detectionResponseFilename = 'PC_spectral_response_CZT0.25x0.25x1.6.mat'
        cfg.scanner.detectorBinThreshold = [20, 30, 40, 60, 80, 100, 160]
        cfg.scanner.detectorSumBins = 0

        cfg.sim.enableQuantumNoise = 1
        cfg.sim.startViewId = 0
        cfg.sim.stopViewId = 2
        cfg.sim.subViewCount = 1

        cfg.det.totalNumCells = 5
        #cfg.det.cosBetas = np.ones([cfg.det.totalNumCells, 1], dtype=np.single)

        cfg = feval(cfg.protocol.spectrumCallback, cfg)

        cfg.thisSubView = np.float32(np.random.random([cfg.det.totalNumCells, cfg.sim.Evec.size]) * 100)
        cfg.thisSubView[1:4, 3:5] = 0

        viewId = 0
        subViewId = 0
        cfg.sim.eNoise = 3500

        assert cfg.sim.Evec is not None
        assert cfg.thisSubView is not None

        cfg = Detection_PC(cfg, viewId, subViewId)

        assert cfg.sim.Wvec is not None
        assert cfg.thisView is not None


if __name__ == "__main__":
    unittest.main()
