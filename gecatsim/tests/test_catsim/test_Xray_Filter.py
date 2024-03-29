from gecatsim.pyfiles.CommonTools import *
import unittest.mock
from gecatsim.pyfiles import CommonTools
from gecatsim.pyfiles.Xray_Filter import Xray_Filter

class Test_Xray_Filter(unittest.TestCase):

    def test_Xray_Filter(self):
        cfg = CommonTools.CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                              "../examples/cfg/Protocol_Sample_axial")
        cfg.sim.subViewCount = 1

        cfg = feval(cfg.scanner.detectorCallback, cfg)
        cfg = feval(cfg.scanner.focalspotCallback, cfg)
        cfg = feval(cfg.physics.rayAngleCallback, cfg)
        cfg = feval(cfg.protocol.spectrumCallback, cfg)

        cfg.protocol.bowtie = 'medium'
        cfg.protocol.flatFilter = ['al', 0.1, 'water', 2]

        cfg = Xray_Filter(cfg)
        trans = cfg.src.filterTrans.reshape(cfg.scanner.detectorColCount, cfg.scanner.detectorRowCount, cfg.spec.nEbin)

        self.assertAlmostEqual(trans.max(), 0.9266918, places=4)
        assert trans.shape == (900, 16, 20)
        assert trans is not None
