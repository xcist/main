import unittest
from unittest.mock import mock_open, patch
from gecatsim.pyfiles.readViewWeighting import readViewWeighting

class TestReadViewWeighting(unittest.TestCase):
    def test_read_view_weighting_file_exists(self):
        mock_file_content = """0.5 0.0 1.0 1.2 0.3 0.825 0.175 0.0 0.0 0.0
                               0.6 0.0 1.5 1.3 0.4 0.825 0.175 0.0 0.0 0.0"""
        with patch('builtins.open', mock_open(read_data=mock_file_content)):
            vw = readViewWeighting(1.0, 1.5)
            self.assertEqual(vw['collimation'], 0.6)
            self.assertEqual(vw['pitch'], 1.5)
            self.assertEqual(vw['zsf'], 1.0)
            self.assertEqual(vw['kw'], 0.0)
            self.assertEqual(vw['beta_0'], 0.825)
            self.assertEqual(vw['beta_t'], 0.175)
            self.assertEqual(vw['vct_k'], 0.0)
            self.assertEqual(vw['vct_q'], 0.0)
            self.assertEqual(vw['vct_r'], 0.0)

    def test_read_view_weighting_file_not_exists(self):
        with patch('builtins.open', side_effect=FileNotFoundError):
            vw = readViewWeighting(1.0, 1.5)
            self.assertEqual(vw['beta_0'], 0.825)
            self.assertEqual(vw['beta_t'], 0.175)
            self.assertEqual(vw['vct_k'], 0.0)
            self.assertEqual(vw['vct_q'], 0.0)
            self.assertEqual(vw['vct_r'], 0.0)
            self.assertEqual(vw['zsf'], 1.0)
            self.assertEqual(vw['kw'], 0.0)

if __name__ == '__main__':
    unittest.main()