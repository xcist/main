import unittest
from unittest.mock import *
import numpy as np
from gecatsim.pyfiles.CommonTools import CFG
from gecatsim.pyfiles.HVLDoseCalculations import HVLDoseCalculations
import unittest
from unittest.mock import MagicMock
import numpy as np
from gecatsim.pyfiles.CommonTools import CFG
from gecatsim.pyfiles.HVLDoseCalculations import HVLDoseCalculations

class TestHVLDoseCalculations(unittest.TestCase):

    def test_HVLDoseCalculations(self):
        cfg = CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                  "../examples/cfg/Protocol_Sample_axial", "../examples/cfg/Scanner_Sample_generic.cfg")

        cfg.dose_distance = 1000
        cfg.col_count = 2
        cfg.row_count = 2
        cfg.calculate_central_dose = True
        cfg.calculate_hvl = True

        spec = MagicMock()
        spec.Evec = np.array([[10, 20, 30], [40, 50, 60]])
        spec.netIvec = np.ones((3, 4))

        HVLDoseCalculations(cfg, spec)


if __name__ == '__main__':
    unittest.main()