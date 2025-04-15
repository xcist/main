import unittest
from unittest.mock import MagicMock
import numpy as np
from gecatsim.pyfiles.CommonTools import CFG
from gecatsim.pyfiles.Resample_Spectrum_Bowtie_FlatFilter import Resample_Spectrum_Bowtie_FlatFilter

class TestResampleSpectrumBowtieFlatFilter(unittest.TestCase):

    def test_resample_spectrum_bowtie_flatfilter(self):
        cfg = CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                  "../examples/cfg/Protocol_Sample_axial", "../examples/cfg/Scanner_Sample_generic.cfg")

        cfg.start_angle = 0
        cfg.views_per_rotation = 360
        cfg.rotation_direction = 1
        cfg.mA_modulation = 0.5

        bowtie = MagicMock()
        bowtie.transVec = np.ones((360, 1))
        bowtie.sinalphas = np.ones(360)
        bowtie.singammas = np.ones(360)

        spec = MagicMock()
        spec.Ivec = np.ones((360, 1))
        spec.sinalphas = np.ones(360)
        spec.singammas = np.ones(360)

        FiltrationTransVec = np.ones((360, 1))

        bowtie, spec, FiltrationTransVec = Resample_Spectrum_Bowtie_FlatFilter(cfg, bowtie, spec, FiltrationTransVec)

        self.assertTrue(np.array_equal(spec.netIvec, spec.Ivec * bowtie.transVec * FiltrationTransVec))

if __name__ == '__main__':
    unittest.main()