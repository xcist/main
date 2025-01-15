import unittest
from gecatsim.pyfiles.get_catsim_spectrum import get_catsim_spectrum
from gecatsim.pyfiles.CommonTools import *

class TestGetCatsimSpectrum(unittest.TestCase):
    def test_get_catsim_spectrum(self):

        cfg = CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                  "../examples/cfg/Protocol_Sample_axial")
        cfg.sim.subViewCount = 1
        cfg.det.totalNumCells = 14400
        cfg.spec.nEbin = 20
        cfg.Ivec = [1, 2, 3]
        cfg.netIvec = [0.9, 1.8, 2.7]
        cfg.Evec = [10, 20, 30]

        # Test spectrum before filters
        cfg = feval(cfg.protocol.spectrumCallback, cfg)
        ss, spec, cfg = get_catsim_spectrum(cfg, sim_step=0)
        self.assertIsNotNone(ss)
        self.assertIsNotNone(spec.Ivec)

        # # Test spectrum after bowtie and flat filter
        # cfg = feval(cfg.scanner.detectorCallback, cfg)
        # cfg = feval(cfg.scanner.focalspotCallback, cfg)
        # cfg = feval(cfg.scanner.bowtieCallback, cfg)
        # cfg = feval(cfg.scanner.filtrationCallback, cfg)
        # ss, spec, cfg = get_catsim_spectrum(cfg, sim_step=1)
        # self.assertIsNotNone(ss)
        # self.assertIsNotNone(spec.netIvec)
        #
        # # Test spectrum after detector prefilter
        # if hasattr(cfg.scanner, 'detectorPrefilterCallback'):
        #     cfg = feval(cfg.scanner.detectorPrefilterCallback, cfg)
        # ss, spec, cfg = get_catsim_spectrum(cfg, sim_step=2)
        # self.assertIsNotNone(ss)
        # self.assertIsNotNone(spec.netIvec)

if __name__ == '__main__':
    unittest.main()