from gecatsim.pyfiles.CommonTools import *
import unittest.mock
from gecatsim.pyfiles import CommonTools
from gecatsim.pyfiles.Spectrum import Spectrum

class Test_Spectrum(unittest.TestCase):
    def test1_Spectrum(self):
        cfg = CommonTools.CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                              "../examples/cfg/Protocol_Sample_axial")

        cfg.det.totalNumCells = 5;

        cfg.protocol.spectrumFilename = "tungsten_tar7.0_120_filt.dat"
        cfg.physics.energyCount = 10
        cfg.protocol.spectrumScaling = 1
        cfg.physics.monochromatic = -1

        cfg.protocol.mA = 200
        cfg.protocol.rotationTime = 1
        cfg.protocol.viewsPerRotation = 984

        cfg.sim.subViewCount = 1
        cfg.protocol.dutyRatio = 1

        cfg = Spectrum(cfg)

        check_value(cfg.spec.nEbin)
        check_value(cfg.spec.Evec)
        check_value(cfg.spec.Ivec)

        assert cfg.spec.nEbin == 10
        assert (cfg.spec.Evec == [  6.,  18.,  30.,  42.,  54.,  66.,  78.,  90., 102., 114.] ).all()
        assert cfg.spec.Ivec is not None



