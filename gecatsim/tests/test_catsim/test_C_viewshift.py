import unittest
from unittest.mock import MagicMock
from gecatsim.pyfiles.C_viewshift import C_viewshift
from gecatsim.pyfiles.CommonTools import *

class TestCViewshift(unittest.TestCase):

    def test_C_viewshift_success(self):
        cfg = CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                  "../examples/cfg/Protocol_Sample_axial")

        cfg.clib.viewshift = MagicMock(return_value=0)
        rows = 2
        cols = 2
        planes = 1
        pmax = 10
        view = 1.0
        ndx = [1, 2]
        coef = [0.5, 0.5]
        output = [0.0, 0.0]

        result = C_viewshift(cfg, rows, cols, planes, pmax, view, ndx, coef, output)
        cfg.clib.viewshift.assert_called_once()
        self.assertEqual(result, output)

    def test_C_viewshift_failure(self):
        cfg = CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                  "../examples/cfg/Protocol_Sample_axial")

        cfg.clib.viewshift = MagicMock(return_value=-1)
        rows = 2
        cols = 2
        planes = 1
        pmax = 10
        view = 1.0
        ndx = [1, 2]
        coef = [0.5, 0.5]
        output = [0.0, 0.0]

        with self.assertRaises(Exception) as context:
            C_viewshift(cfg, rows, cols, planes, pmax, view, ndx, coef, output)
        self.assertTrue("Error occurred while calling C function viewshift." in str(context.exception))

if __name__ == '__main__':
    unittest.main()