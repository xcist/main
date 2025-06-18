import unittest
import numpy as np
import io
from contextlib import redirect_stdout
from gecatsim.pyfiles.CommonTools import CFG
from gecatsim.pyfiles.HVLDoseCalculations import HVLDoseCalculations

class TestHVLDoseCalculations(unittest.TestCase):

    def test_hvl_dose_calculations(self):
        cfg = CFG("../examples/cfg/Phantom_Sample",
                  "../examples/cfg/Scanner_Sample_generic",
                  "../examples/cfg/Protocol_Sample_axial")

        cfg.col_count = 2
        cfg.row_count = 2
        cfg.dose_distance = 1000
        cfg.calculate_central_dose = True
        cfg.calculate_hvl = True

        spec = type('', (), {})()
        spec.Evec = np.linspace(10, 100, 10).reshape(-1, 1)
        spec.netIvec = np.ones((10, 5)) * 1e5  # mock photon counts

        # Capture printed output
        f = io.StringIO()
        with redirect_stdout(f):
            HVLDoseCalculations(cfg, spec)
        output = f.getvalue()

        # Assertions based on expected output patterns
        self.assertIn("bCT Air kerma at", output)
        self.assertIn("HVL = ", output)
        self.assertIn("1/4 HVL = ", output)
        self.assertIn("1/8 HVL = ", output)
        self.assertIn("1/10 HVL = ", output)
        self.assertIn("Effective air mu for spectra", output)

        lines = output.splitlines()
        hvl_values = [float(line.split('=')[1].split()[0]) for line in lines if 'HVL =' in line or '1/' in line]
        self.assertTrue(all(hvl > 0 for hvl in hvl_values), "All HVL values should be positive")
        self.assertTrue(hvl_values == sorted(hvl_values), "HVL values should increase with attenuation level")

if __name__ == "__main__":
    unittest.main()
